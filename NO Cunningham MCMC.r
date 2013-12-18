library(foreign)
library(sandwich)
library(ggplot2)
library(lmtest)
library(lme4)
library(boot)
library(MASS)
library(coefplot)
library(ROCR)
library(stargazer)
library(RColorBrewer)
library(bootstrap)
library(rjags)

setwd("C:/Users/Nathaniel Olin/Dropbox/Schoolwork/PS818/Replication Project")
rm(list=ls())

#############
# FUNCTIONS #
#############

# robust clustered standard errors
robclust <- function(model, data, cluster){
  n <- nrow(data)
  cluster.list <- factor(data[,cluster])
  m <- length(levels(cluster.list))
  k <- length(coef(model))
  dfc <- (m/(m-1))*((n-1)/(n-k))
  u <- estfun(model)
  u.clust <- matrix(NA, nrow=m, ncol=k)
  for(j in 1:k){u.clust[,j] <- tapply(u[,j], cluster.list, sum)}
  rob.clust <- vcov(model)%*%(dfc*t(u.clust)%*%u.clust)%*%vcov(model)
}

########
# DATA #
########

# Import replication data
RepData <- read.dta("AJPS_Cunningham Rep data.dta")

# Rename a few problem variables
names(RepData)[c(16:18,22:24)] <- 
  c("spline1_cwonset","spline2_cwonset","spline3_cwonset",
    "spline1_cw","spline2_cw","spline3_cw")

# Create new dataframes and remove missing data (necessary for SEs)
OnsetData <- subset(RepData, ongoingcivilwar!=1)
attach(OnsetData)
OnsetClean <- na.omit(data.frame(kgcid,country,group,year,civilwaronset,
              logfactions,prevconcessions_l,democracy,kin,yrnocwonset,
			  spline1_cwonset,spline2_cwonset,spline3_cwonset))
detach(OnsetData)
attach(RepData)
IncidenceClean <- na.omit(data.frame(kgcid,country,group,year,
                  acdcivilwar1,logfactions,prevconcessions_l,democracy,
                  kin,yrsnocivilwar,spline1_cw,spline2_cw,spline3_cw))
detach(RepData)

#########################
# REPLICATE BASE MODELS #
#########################

# logit model of civil war onset
logOnset <- glm(
  civilwaronset ~ logfactions + prevconcessions_l + democracy + kin + 
  yrnocwonset + spline1_cwonset + spline2_cwonset + spline3_cwonset,
  data=OnsetClean,family=binomial(link="logit"))
logOnsetClust <- robclust(logOnset, OnsetClean, "kgcid")
logOnsetTest <- coeftest(logOnset,logOnsetClust)

# logit model of civil war incidence
logIncidence <- glm(
  acdcivilwar1 ~ logfactions + prevconcessions_l + democracy + kin + 
  yrsnocivilwar + spline1_cw + spline2_cw + spline3_cw,
  data=IncidenceClean,family=binomial(link="logit"))
logIncidenceClust <- robclust(logIncidence, IncidenceClean, "kgcid")
logIncidenceTest <- coeftest(logIncidence,logIncidenceClust)

###############
# SPATIAL LAG #
###############

library(ape)
# Connectivity matrix for incidence

N <- nrow(OnsetClean)
W1 <- matrix(0,nrow=N,ncol=N)
for(i in 1:N){
  tc <- which(OnsetClean$country==OnsetClean$country[i])
  tyl <- which(OnsetClean$year>=OnsetClean$year[i]-5)
  tyh <- which(OnsetClean$year<OnsetClean$year[i])
  ty <- intersect(tyl,tyh)
  W1[i,intersect(tc,ty)] <- 1/(length(intersect(tc,ty)))
  W1[i,i] <- 0
}
OnsetClean$spacelag.onset <- as.vector(W1%*%OnsetClean$civilwaronset)
OnsetAutoTest <- Moran.I(OnsetClean$civilwaronset,W1)
#ggplot(OnsetClean) + 
#  geom_density(aes(x=spacelag.onset),fill="light blue") +
#  theme_minimal() + 
#  labs(x="Proportion of neighbors going to war",y="Density")

N <- nrow(IncidenceClean)
W2 <- matrix(0,nrow=N,ncol=N)
for(i in 1:N){
  tc <- which(IncidenceClean$country==IncidenceClean$country[i])
  ty <- which(IncidenceClean$year==IncidenceClean$year[i]-1)
  W2[i,intersect(tc,ty)] <- 1/length(intersect(tc,ty))
  W2[i,i] <- 0
}
IncidenceClean$spacelag.incidence <- 
  as.vector(W2%*%IncidenceClean$acdcivilwar1)
IncidenceAutoTest <- Moran.I(IncidenceClean$acdcivilwar1,W2)

MoranTable <- 
 rbind(rbind(OnsetAutoTest[1:4]),rbind(IncidenceAutoTest[1:4]))
rownames(MoranTable) <- c("Onset","Incidence")
colnames(MoranTable) <- c("Observed","Expected","SD","p value") 
sink("MoranTable.tex")
  print(xtable(MoranTable,digits=4),floating=FALSE,
    math.style.negative=TRUE,booktabs=TRUE)
sink()

IncidenceClean$incplot <- as.factor(IncidenceClean$acdcivilwar1)
levels(IncidenceClean$incplot) <- 
  c("Observations at peace","Observations at war")

mycolors <- brewer.pal(2,"Set1")[c(2,1)]
pdf("densityIncidence.pdf")
ggplot(IncidenceClean,aes(x=spacelag.incidence,fill=incplot,color=incplot)) + 
  geom_density(alpha=0.5) +
  scale_fill_manual(values=mycolors,name="") + 
  scale_color_manual(values=mycolors,name="") + 
  #scale_fill_brewer(palette="Set1",name="") + 
  theme_minimal() +
  labs(x="Proportion of neighbors at war",y="Density") +
  theme(legend.position="bottom") +
  ggtitle("Proportion of neighboring dyads in conflict")
dev.off()

#mycolors <- brewer.pal(4,"Set1")[c(2,4)]
#ggplot(IncidenceClean,aes(x=spacelag.incidence,color=as.factor(acdcivilwar1))) + 
#  geom_density(size=1) +
#  scale_color_manual(values=mycolors) +
#  labs(x="Proportion of neighbors at war",y="Density") +
#  theme_minimal()


############################
# MARKOV CHAIN MONTE-CARLO #		 
############################

# Onset model and priors
OnsetModel <- "model {
  for(i in 1:length(civilwaronset)) {
    civilwaronset[i] ~ dbern(p[i])
    p[i] <- 1/(1+exp(-z[i]))
    z[i] <- 
	  alpha + beta1*logfactions[i] + beta2*prevconcessions_l[i] + beta3*democracy[i] + 
	  beta4*kin[i] + beta5*yrnocwonset[i] + beta6*spline1_cwonset[i] + 
	  beta7*spline2_cwonset[i] + beta8*spline3_cwonset[i] + gamma*spacelag.onset[i]
  }
  alpha ~ dnorm(0,0.0001)
  beta1 ~ dnorm(0,0.0001)
  beta2 ~ dnorm(0,0.0001)
  beta3 ~ dnorm(0,0.0001)
  beta4 ~ dnorm(0,0.0001)
  beta5 ~ dnorm(0,0.0001) 
  beta6 ~ dnorm(0,0.0001)
  beta7 ~ dnorm(0,0.0001)
  beta8 ~ dnorm(0,0.0001)
  gamma ~ dnorm(0,0.0001)
}"
# Data for sampler
attach(OnsetClean)
OnsetList <- list(
  "civilwaronset"=civilwaronset,"logfactions"=logfactions,
  "prevconcessions_l"=prevconcessions_l,"democracy"=democracy,"kin"=kin,
  "yrnocwonset"=yrnocwonset,"spline1_cwonset"=spline1_cwonset,
  "spline2_cwonset"=spline2_cwonset,"spline3_cwonset"=spline3_cwonset,
  "spacelag.onset"=spacelag.onset)
detach(OnsetClean)

parameters <- c(
  "alpha","beta1","beta2","beta3","beta4","beta5","beta6","beta7","beta8","gamma")
				
# Running the model
OnsetMCMC <- jags.model(
  textConnection(OnsetModel),data=OnsetList,n.chains = 4, n.adapt=500)
update(OnsetMCMC,20000); # Burnin
OnsetSamples <- coda.samples(
  OnsetMCMC,variable.names=parameters,n.iter=40000)
OnsetMCMC.out <- summary(OnsetSamples)

# Table output
OnsetNames <- c(
  "(Intercept)","Log SD factions","Previous concessions","Democracy",
  "Kin","Years since civil war onset","Spline 1","Spline 2","Spline 3",
  "Neighbor effect")

OnsetTable <- round(cbind(
  c(logOnsetTest[,"Estimate"],NA),c(logOnsetTest[,"Std. Error"],NA),
  OnsetMCMC.out[[1]][,"Mean"],OnsetMCMC.out[[1]][,"SD"]),2)
rownames(OnsetTable) <- OnsetNames
colnames(OnsetTable) <- c("Mean","SD","Mean","SD")

sink("OnsetSpace.tex")
  print(xtable(OnsetTable),floating=FALSE,math.style.negative=TRUE,
    booktabs=TRUE)
sink()
  
# Incidence model and priors
IncidenceModel <- "model {
  for (i in 1:length(acdcivilwar1)) {
    acdcivilwar1[i] ~ dbern(p[i])
    p[i] <- 1/(1+exp(-z[i]))
    z[i] <- 
	  alpha + beta1*logfactions[i] + beta2*prevconcessions_l[i] + beta3*democracy[i] + beta4*kin[i] + beta5*yrsnocivilwar[i] + beta6*spline1_cw[i] + beta7*spline2_cw[i] + 
	  beta8*spline3_cw[i] + gamma*spacelag.incidence[i]
  }
  alpha ~ dnorm(0,0.0001)
  beta1 ~ dnorm(0,0.0001)
  beta2 ~ dnorm(0,0.0001)
  beta3 ~ dnorm(0,0.0001)
  beta4 ~ dnorm(0,0.0001)
  beta5 ~ dnorm(0,0.0001) 
  beta6 ~ dnorm(0,0.0001)
  beta7 ~ dnorm(0,0.0001)
  beta8 ~ dnorm(0,0.0001)
  gamma ~ dnorm(0,0.0001)
}"

attach(IncidenceClean)
IncidenceList <- list(
  "acdcivilwar1"=acdcivilwar1,"logfactions"=logfactions,
  "prevconcessions_l"=prevconcessions_l,"democracy"=democracy,"kin"=kin,
  "yrsnocivilwar"=yrsnocivilwar,"spline1_cw"=spline1_cw,"spline2_cw"=spline2_cw,
  "spline3_cw"=spline3_cw,"spacelag.incidence"=spacelag.incidence)
detach(IncidenceClean)

parameters <- c(
  "alpha","beta1","beta2","beta3","beta4","beta5","beta6","beta7","beta8","gamma")
				
# Running the model
IncidenceMCMC <- jags.model(
  textConnection(IncidenceModel),data=IncidenceList,n.chains = 4, n.adapt=500)
update(IncidenceMCMC,20000); # Burnin
IncidenceSamples <- coda.samples(
  IncidenceMCMC,variable.names=parameters,n.iter=40000)

IncidenceMCMC.out <- summary(IncidenceSamples)

# Table output
IncidenceNames <- c(
  "(Intercept)","Log SD factions","Previous concessions","Democracy",
  "Kin","Years since civil war","Spline 1","Spline 2","Spline 3",
  "Neighbor effect")

IncidenceTable <- round(cbind(
  c(logIncidenceTest[,"Estimate"],NA),
  c(logIncidenceTest[,"Std. Error"],NA),
  IncidenceMCMC.out[[1]][,"Mean"],IncidenceMCMC.out[[1]][,"SD"]),2)
rownames(IncidenceTable) <- IncidenceNames
colnames(IncidenceTable) <- c("Mean","SD","Mean","SD")

sink("IncidenceSpace.tex")
  print(xtable(IncidenceTable),floating=FALSE,math.style.negative=TRUE,
    booktabs=TRUE)
sink()

################
# PLOT RESULTS #
################

# Sequence of factions
factions.seq <- seq(1,20,0.1)

# ONSET

# Set up design matrices
attach(OnsetClean)
# Base model
OnsetX <- cbind(1,log(factions.seq),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrnocwonset),mean(spline1_cwonset),
  mean(spline2_cwonset),mean(spline3_cwonset))
B.Onset <- coef(logOnset)
predprob.Onset <- as.vector(plogis(B.Onset%*%t(OnsetX)))
# Spacelag
OnsetSpaceX <- cbind(
  1,log(factions.seq),median(prevconcessions_l),median(democracy),
  median(kin),mean(yrnocwonset),mean(spline1_cwonset),
  mean(spline2_cwonset),mean(spline3_cwonset),mean(spacelag.onset))
B.OnsetSpace <- OnsetMCMC.out[[1]][,"Mean"]
predprob.OnsetSpace <- as.vector(plogis(B.OnsetSpace%*%t(OnsetSpaceX)))
detach(OnsetClean)

probs.Onset <- data.frame(cbind(
  factions.seq,predprob.Onset,predprob.OnsetSpace))
names(probs.Onset) <- c("factions.seq","base","space")

OnsetPlot <- 
  ggplot(probs.Onset) + 
  geom_line(aes(factions.seq,base,color="Base model"),size=1) +
  geom_line(aes(factions.seq,space,color="Spatial lag"),size=1)
pdf("probsOnsetSpace.pdf",width=5,height=5)
OnsetPlot + theme_minimal() + 
  scale_y_continuous(limits=c(0,0.6)) +
  scale_color_brewer(palette="Set1",name="Model") +
  labs(x="Number of factions",y="") +
  theme(legend.position=c(1,0),legend.justification=c(1,0),
        legend.background = element_rect()) +
  ggtitle("Predicted probability of civil war onset")
dev.off()


# INCIDENCE

# Set up design matrices
attach(IncidenceClean)
# Base model
IncidenceX <- cbind(1,log(factions.seq),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrsnocivilwar),mean(spline1_cw),
  mean(spline2_cw),mean(spline3_cw))
B.Incidence <- coef(logIncidence)
predprob.Incidence<- as.vector(plogis(B.Incidence%*%t(IncidenceX)))
# Spacelag
IncidenceSpaceX <- cbind(
  1,log(factions.seq),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrsnocivilwar),mean(spline1_cw),
  mean(spline2_cw),mean(spline3_cw),mean(spacelag.incidence))
#B.IncidenceSpace <- coef(logIncidence6)
B.IncidenceSpace <- IncidenceMCMC.table[,1]
predprob.IncidenceSpace <- as.vector(plogis(B.IncidenceSpace%*%t(IncidenceSpaceX)))
detach(IncidenceClean)

probs.Incidence <- data.frame(cbind(
  factions.seq,predprob.Incidence,predprob.IncidenceSpace))
names(probs.Incidence) <- c("factions.seq","base","space")

IncidencePlot <- 
  ggplot(probs.Incidence) + 
  geom_line(aes(factions.seq,base,color="Base model"),size=1) +
  geom_line(aes(factions.seq,space,color="Spatial lag"),size=1)
pdf("probsIncidenceSpace.pdf",width=5,height=5)
IncidencePlot + theme_minimal() + 
  scale_y_continuous(limits=c(0,0.6)) +
  scale_color_brewer(palette="Set1",name="Model") +
  labs(x="Number of factions",y="") +
  theme(legend.position=c(1,0),legend.justification=c(1,0),
        legend.background = element_rect()) +
  ggtitle("Predicted probability of civil war incidence")
dev.off()

############
# REAL TRY #
############

space.seq <- seq(0,1,0.01)
sample.seq <- seq(1,40000,10)

attach(OnsetClean)
OnsetX.2 <- cbind(1,log(2),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrnocwonset),mean(spline1_cwonset),
  mean(spline2_cwonset),mean(spline3_cwonset),space.seq)
detach(OnsetClean)
B.Onset <- OnsetSamples[[1]][sample.seq,]
#B.Onset <- mvrnorm(10000,coef(logOnset6),logOnset6Clust)
predprob.Onset <- plogis(B.Onset%*%t(OnsetX.2))
predprob.Onset <- t(apply(predprob.Onset,2,quantile,c(0.025,0.925)))
predprob.Onset <- data.frame(cbind(
  space.seq,predprob.Onset,apply(predprob.Onset,1,mean)))
colnames(predprob.Onset) <- c("factions.seq","ci.low","ci.high","mean")

attach(IncidenceClean)
IncidenceX.2 <- cbind(1,log(2),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrsnocivilwar),mean(spline1_cw),
  mean(spline2_cw),mean(spline3_cw),space.seq)
detach(IncidenceClean)
B.Incidence <- IncidenceSamples[[1]][sample.seq,]
#B.Incidence <- mvrnorm(10000, coef(logIncidence), logIncidenceClust)
predprob.Incidence <- plogis(B.Incidence%*%t(IncidenceX.2))
predprob.Incidence <- 
  t(apply(predprob.Incidence,2,quantile,c(0.025,0.925)))
predprob.Incidence <- data.frame(cbind(
  space.seq,predprob.Incidence,apply(predprob.Incidence,1,mean)))
colnames(predprob.Incidence) <- c("factions.seq","ci.low","ci.high","mean")

# Plot simulated confidence intervals
OnsetPlot <- 
  ggplot(predprob.Onset, aes(space.seq,mean)) + geom_line() + 
  geom_ribbon(aes(ymin=ci.low,ymax=ci.high),alpha=0.3,fill="blue")
IncidencePlot <- 
  ggplot(predprob.Incidence, aes(space.seq,mean)) + geom_line() + 
  geom_ribbon(aes(ymin=ci.low,ymax=ci.high),alpha=0.3,fill="blue")

pdf("finalspaceOnset2.pdf",width=5,height=5)
OnsetPlot + theme_minimal() + 
  scale_x_continuous(limits=c(0,1),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,0.6)) +
  labs(x="Proportion of neighboring dyads at war",y="") +
  ggtitle("Predicted probability of civil war onset")
dev.off()
pdf("finalspaceIncidence.pdf",width=5,height=5)
IncidencePlot + theme_minimal() + 
  scale_x_continuous(limits=c(0,1),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,0.6)) +
  labs(x="Proportion of neighboring dyads at war",y="") +
  ggtitle("Predicted probability of civil war incidence")
dev.off()





########################################3













# ROC CURVES

predOnset <- 
  prediction(logOnset$fitted.values,OnsetClean$civilwaronset)
performOns <- performance(predOnset,"tpr","fpr")
tpr <- unlist(performOns@y.values)
fpr <- unlist(performOns@x.values)
rocOns <- data.frame(tpr,fpr)

predOnsetSpace <- 
  prediction(logOnset6$fitted.values,OnsetClean$civilwaronset)
performOnsSpace <- performance(predOnsetSpace,"tpr","fpr")
tprSpace <- unlist(performOnsSpace@y.values)
fprSpace <- unlist(performOnsSpace@x.values)
rocOnsSpace <- data.frame(tprSpace,fprSpace)

pdf("rocOnsSpace.pdf",width=5,height=5)
ggplot(rocOns) + 
  geom_line(aes(x=fpr,y=tpr,color="Base model")) + 
  geom_line(data=rocOnsSpace,aes(x=fprSpace,y=tprSpace,
    color="Spatial lag")) +
  geom_abline(intercept=0,slope=1,linetype="dashed") +
  theme_minimal() + 
  theme(legend.position=c(1,0),legend.justification=c(1,0),
    legend.background=element_rect()) +
  scale_color_brewer(palette="Set1",name="Models") +
  labs(x="False positive rate",y="True positive rate") +
  ggtitle("Base onset model vs spatial lag model")
dev.off()


predIncidence <- 
  prediction(logIncidence$fitted.values,IncidenceClean$acdcivilwar1)
performInc <- performance(predIncidence,"tpr","fpr")
tpr <- unlist(performInc@y.values)
fpr <- unlist(performInc@x.values)
rocInc <- data.frame(tpr,fpr)

predIncidenceSpace <- 
  prediction(logIncidence6$fitted.values,IncidenceClean$acdcivilwar1)
performIncSpace <- performance(predIncidenceSpace,"tpr","fpr")
tprSpace <- unlist(performIncSpace@y.values)
fprSpace <- unlist(performIncSpace@x.values)
rocIncSpace <- data.frame(tprSpace,fprSpace)


pdf("rocIncSpace.pdf",width=5,height=5)
ggplot(rocInc) + 
  geom_line(aes(x=fpr,y=tpr,color="Base model")) + 
  geom_line(data=rocIncSpace,aes(x=fprSpace,y=tprSpace,
    color="Spatial lag")) +
  geom_abline(intercept=0,slope=1,linetype="dashed") +
  theme_minimal() + 
  theme(legend.position=c(1,0),legend.justification=c(1,0),
    legend.background=element_rect()) +
  scale_color_brewer(palette="Set1",name="Models") +
  labs(x="False positive rate",y="True positive rate") +
  ggtitle("Base incidence model vs spatial lag model")
dev.off()
