library(foreign)
library(sandwich)
library(ggplot2)
library(lmtest)
library(lme4)
library(boot)
library(MASS)
library(coefplot)
library(ROCR)
library(ggROC)
library(stargazer)
library(RColorBrewer)
library(bootstrap)

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
names(RepData)
attributes(RepData)$var.labels

# Rename a few problem variables
names(RepData)[c(16:18,22:24)] <- c(
  "spline1_cwonset","spline2_cwonset","spline3_cwonset", "spline1_cw",
  "spline2_cw","spline3_cw")

# Create new dataframes and remove missing data (necessary for SEs)
OnsetData <- subset(RepData, ongoingcivilwar!=1)
attach(OnsetData)
OnsetClean <- na.omit(data.frame(
  kgcid,country,group,year,civilwaronset,logfactions,prevconcessions_l,
  democracy,kin,yrnocwonset,spline1_cwonset,spline2_cwonset,
  spline3_cwonset))
detach(OnsetData)
attach(RepData)
IncidenceClean <- na.omit(data.frame(
  kgcid,country,group,year,acdcivilwar1,logfactions,prevconcessions_l,democracy,kin,yrsnocivilwar,spline1_cw,spline2_cw,spline3_cw))
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

#########################
# NEW ROBUSTNESS CHECKS #
#########################

# TEST WITH LAGGED # SD MOVEMENTS IN PREVIOUS YEAR

# Onset: define variable
N <- nrow(OnsetClean)
OnsetClean$lagmovements <- rep(0,N)
for(i in 1:N){
  # Identify observations in the same country
  tc <- which(OnsetClean$country==OnsetClean$country[i])
  # Identify observations in the previous year
  ty <- which(OnsetClean$year==OnsetClean$year[i]-1)
  # Identify where these overlap and save the count
  OnsetClean$lagmovements[i] <- length(intersect(tc,ty))
}
# Test onset with lagged number of SD movements
logOnset2 <- glm(
  civilwaronset ~ logfactions + prevconcessions_l + democracy + kin + 
  yrnocwonset + spline1_cwonset + spline2_cwonset + spline3_cwonset + 
  lagmovements, data=OnsetClean, family=binomial(link="logit"))
logOnset2Clust <- robclust(logOnset2, OnsetClean, "kgcid")
logOnset2Test <- coeftest(logOnset2,logOnset2Clust)

# Incidence: define variable
N <- nrow(IncidenceClean)
IncidenceClean$lagmovements <- rep(0,N)
for(i in 1:N){
  # Identify observations in the same country
  tc <- which(IncidenceClean$country==IncidenceClean$country[i])
  # Identify observations in the previous year
  ty <- which(IncidenceClean$year==IncidenceClean$year[i]-1)
  # Identify where these overlap and save the count
  IncidenceClean$lagmovements[i] <- length(intersect(tc,ty))
}
# Test incidence with lagged number of SD factions
logIncidence2 <- glm(
  acdcivilwar1 ~ logfactions + prevconcessions_l + democracy + kin + 
  yrsnocivilwar + spline1_cw + spline2_cw + spline3_cw + lagmovements,
  data=IncidenceClean,family=binomial(link="logit"))
logIncidence2Clust <- robclust(logIncidence2,IncidenceClean,"kgcid")
logIncidence2Test <- coeftest(logIncidence2,logIncidence2Clust)

# COUNTRY DUMMIES?
#logOnset3 <- glm(
#  civilwaronset ~ logfactions + prevconcessions_l + democracy + kin + 
#  yrnocwonset + spline1_cwonset + spline2_cwonset + spline3_cwonset + 
#  as.factor(country),data=OnsetClean,family=binomial(link="logit"))
#logOnset3Clust <- robclust(logOnset3, OnsetClean, "kgcid")
#logOnset3Test <- coeftest(logOnset3,logOnset3Clust)

#logIncidence3 <- glm(
#  acdcivilwar1 ~ logfactions + prevconcessions_l + democracy + kin + yrsnocivilwar + spline1_cw + spline2_cw + spline3_cw + as.factor(country),data=IncidenceClean,family=binomial(link="logit"))
#logIncidence3Clust <- robclust(logIncidence3, IncidenceClean, "kgcid")
#logIncidence3Test <- coeftest(logIncidence3,logIncidence3Clust)

# PAST CIVIL WAR IN COUNTRY

OnsetClean$pastcwonset <- rep(0,nrow(OnsetClean))
for(i in 1:nrow(OnsetClean)){
  ty <- which(OnsetClean$year<OnsetClean$year[i])
  tc <- which(OnsetClean$country==OnsetClean$country[i])
  past <- nrow(subset(OnsetClean[intersect(ty,tc),],civilwaronset==1))
  OnsetClean$pastcwonset[i] <- ifelse(past>0, 1, 0)
}

logOnset3 <- glm(
  civilwaronset ~ logfactions + prevconcessions_l + democracy + kin + 
  yrnocwonset + spline1_cwonset + spline2_cwonset + spline3_cwonset + 
  pastcwonset,data=OnsetClean,family=binomial(link="logit"))
logOnset3Clust <- robclust(logOnset3, OnsetClean, "kgcid")
logOnset3Test <- coeftest(logOnset3,logOnset3Clust)

IncidenceClean$pastcw <- rep(0,nrow(IncidenceClean))
for(i in 1:nrow(IncidenceClean)){
  ty <- which(IncidenceClean$year<IncidenceClean$year[i])
  tc <- which(IncidenceClean$country==IncidenceClean$country[i])
  past <- nrow(subset(IncidenceClean[intersect(ty,tc),],acdcivilwar1==1))
  IncidenceClean$pastcw[i] <- ifelse(past>0, 1, 0)
}

logIncidence3 <- glm(
  acdcivilwar1 ~ logfactions + prevconcessions_l + democracy + kin + yrsnocivilwar + spline1_cw + spline2_cw + spline3_cw + pastcw,data=IncidenceClean,family=binomial(link="logit"))
logIncidence3Clust <- robclust(logIncidence3, IncidenceClean, "kgcid")
logIncidence3Test <- coeftest(logIncidence3,logIncidence3Clust)


# NO, RANDOM EFFECTS!
logOnset4 <- glmer(
  civilwaronset ~ logfactions + prevconcessions_l + democracy + kin + 
  yrnocwonset + spline1_cwonset + spline2_cwonset + spline3_cwonset + 
  (1|country),data=OnsetClean,family=binomial(link="logit"))
logIncidence4 <- glmer(
  acdcivilwar1 ~ logfactions + prevconcessions_l + democracy + kin + 
  yrsnocivilwar + spline1_cw + spline2_cw + spline3_cw + (1 | country),
  data=IncidenceClean,family=binomial(link="logit"))

# KITCHEN SINK

logOnset5 <- glmer(
  civilwaronset ~ logfactions + prevconcessions_l + democracy + kin + 
  yrnocwonset + spline1_cwonset + spline2_cwonset + spline3_cwonset + 
  lagmovements + pastcwonset + (1|country),
  data=OnsetClean,family=binomial(link="logit"))
logIncidence5 <- glmer(
  acdcivilwar1 ~ logfactions + prevconcessions_l + democracy + kin + 
  yrsnocivilwar + spline1_cw + spline2_cw + spline3_cw + lagmovements +
  + pastcw + (1 | country),
  data=IncidenceClean,family=binomial(link="logit"))

##########################
# ROBUSTNESS CHECK TABLE #
##########################

OnsetNames <- c(
  "Log SD factions","Previous concessions","Democracy","Kin",
  "Years since civil war onset","Spline 1","Spline 2","Spline 3",
  "Lagged SD movements","Past civil war")
IncidenceNames <- c(
  "Log SD factions","Previous concessions","Democracy","Kin",
  "Years since civil war","Spline 1","Spline 2","Spline 3",
  "Lagged SD movements","Past civil war")
OnsetSE <- list(
  logOnsetTest[,"Std. Error"],
  logOnset2Test[,"Std. Error"],
  logOnset3Test[,"Std. Error"])
OnsetP <- list(
  logOnsetTest[,"Pr(>|z|)"],
  logOnset2Test[,"Pr(>|z|)"],
  logOnset3Test[,"Pr(>|z|)"])
IncidenceSE <- list(
  logIncidenceTest[,"Std. Error"],
  logIncidence2Test[,"Std. Error"],
  logIncidence3Test[,"Std. Error"])
IncidenceP <- list(
  logIncidenceTest[,"Pr(>|z|)"],
  logIncidence2Test[,"Pr(>|z|)"],
  logIncidence3Test[,"Pr(>|z|)"])

# Big tables of numbers
sink("OnsetTable.tex")
stargazer(logOnset,logOnset2,logOnset3,logOnset4,logOnset5,
  se=OnsetSE,p=OnsetP,dep.var.labels="Civil war onset",
  covariate.labels=OnsetNames,model.names=FALSE,
  omit=11:77,omit.stat="bic",digits=2,align=TRUE,
  float=FALSE,header=FALSE,style="ajps")
sink()
sink("IncidenceTable.tex")
stargazer(logIncidence,logIncidence2,logIncidence3,logIncidence4,
  logIncidence5,se=IncidenceSE,p=IncidenceP,omit=11:77,omit.stat="bic",
  dep.var.labels="Civil war incidence",covariate.labels=IncidenceNames,
  model.names=FALSE,digits=2,align=TRUE,float=FALSE,style="ajps")
sink()

##############
# ROC CURVES #
##############

# Civil war onset
predOnset <- 
  prediction(logOnset$fitted.values,OnsetClean$civilwaronset)
performOns <- performance(predOnset,"tpr","fpr")
tpr <- unlist(performOns@y.values)
fpr <- unlist(performOns@x.values)
rocOns <- data.frame(tpr,fpr)

predOnset2 <- 
  prediction(logOnset2$fitted.values,OnsetClean$civilwaronset)
performOns2 <- performance(predOnset2,"tpr","fpr")
tpr2 <- unlist(performOns2@y.values)
fpr2 <- unlist(performOns2@x.values)
rocOns2 <- data.frame(tpr2,fpr2)

predOnset3 <- 
  prediction(logOnset3$fitted.values,OnsetClean$civilwaronset)
performOns3 <- performance(predOnset3,"tpr","fpr")
tpr3 <- unlist(performOns3@y.values)
fpr3 <- unlist(performOns3@x.values)
rocOns3 <- data.frame(tpr3,fpr3)

predOnset4 <- 
  prediction(fitted(logOnset4),OnsetClean$civilwaronset)
performOns4 <- performance(predOnset4,"tpr","fpr")
tpr4 <- unlist(performOns4@y.values)
fpr4 <- unlist(performOns4@x.values)
rocOns4 <- data.frame(tpr4,fpr4)

predOnset5 <- 
  prediction(fitted(logOnset5),OnsetClean$civilwaronset)
performOns5 <- performance(predOnset5,"tpr","fpr")
tpr5 <- unlist(performOns5@y.values)
fpr5 <- unlist(performOns5@x.values)
rocOns5 <- data.frame(tpr5,fpr5)

# Civil war incidence
predIncidence <- 
  prediction(logIncidence$fitted.values,IncidenceClean$acdcivilwar1)
performInc <- performance(predIncidence, "tpr","fpr")
tpr <- unlist(performInc@y.values)
fpr <- unlist(performInc@x.values)
rocInc <- data.frame(tpr,fpr)

predIncidence2 <- 
  prediction(logIncidence2$fitted.values,IncidenceClean$acdcivilwar1)
performInc2 <- performance(predIncidence2, "tpr","fpr")
tpr2 <- unlist(performInc2@y.values)
fpr2 <- unlist(performInc2@x.values)
rocInc2 <- data.frame(tpr2,fpr2)

predIncidence3 <- 
  prediction(logIncidence3$fitted.values,IncidenceClean$acdcivilwar1)
performInc3 <- performance(predIncidence3, "tpr","fpr")
tpr3 <- unlist(performInc3@y.values)
fpr3 <- unlist(performInc3@x.values)
rocInc3 <- data.frame(tpr3,fpr3)

predIncidence4 <- 
  prediction(fitted(logIncidence4),IncidenceClean$acdcivilwar1)
performInc4 <- performance(predIncidence4, "tpr","fpr")
tpr4 <- unlist(performInc4@y.values)
fpr4 <- unlist(performInc4@x.values)
rocInc4 <- data.frame(tpr4,fpr4)

predIncidence5 <- 
  prediction(fitted(logIncidence5),IncidenceClean$acdcivilwar1)
performInc5 <- performance(predIncidence5, "tpr","fpr")
tpr5 <- unlist(performInc5@y.values)
fpr5 <- unlist(performInc5@x.values)
rocInc5 <- data.frame(tpr5,fpr5)

# Plots
rocsize <- 0

pdf("rocOns.pdf",width=5,height=5)
ggplot(rocOns) + 
  geom_line(aes(x=fpr,y=tpr,color="Base model"),size=rocsize) + 
  geom_line(data=rocOns2,aes(x=fpr2,y=tpr2,
    color="SD movements"),size=rocsize) +
  geom_line(data=rocOns4,aes(x=fpr4,y=tpr4,color="Random effects"),
    size=rocsize) +
  geom_line(data=rocOns3,aes(x=fpr3,y=tpr3,color="Past civil war"),
    size=rocsize) +
  geom_line(data=rocOns5,aes(x=fpr5,y=tpr5,color="With all checks"),
    size=rocsize) +
  geom_line(aes(x=c(0,1),y=c(0,1)),linetype="dashed") +
  theme_minimal() + 
  theme(legend.position=c(1,0),legend.justification=c(1,0),
    legend.background=element_rect()) +
  scale_color_brewer(palette="Set1",name="Models") +
  labs(x="False positive rate",y="True positive rate") +
  ggtitle("Onset model and robustness checks")
dev.off()
  
pdf("rocInc.pdf",width=5,height=5)
ggplot(rocInc) + 
  geom_line(aes(x=fpr,y=tpr,color="Base model"),
    size=rocsize) + 
  geom_line(data=rocInc2,aes(x=fpr2,y=tpr2,color="SD movements"),
    size=rocsize) +
  geom_line(data=rocInc4,aes(x=fpr4,y=tpr4,color="Random effects"),
    size=rocsize) +
  geom_line(data=rocInc3,aes(x=fpr3,y=tpr3,color="Past civil war"),
    size=rocsize) +
  geom_line(data=rocInc5,aes(x=fpr5,y=tpr5,color="With all checks"),
    size=rocsize) +
  geom_line(aes(x=c(0,1),y=c(0,1)),linetype="dashed") +
  theme_minimal() + 
  theme(legend.position=c(1,0),legend.justification=c(1,0),
    legend.background=element_rect()) +
  scale_color_brewer(palette="Set1",name="Models") +
  labs(x="False positive rate",y="True positive rate") +
  ggtitle("Incidence model and robustness checks")
dev.off()

####################
# COEFFICIENT PLOT #
####################

# Onset plot
varnames <- c(
  "Log number of\nself-determination\nfactions",
  "Previous concessions","Democracy","Kin in neighboring\nstate",
  "Years since civil\nwar onset",
  "Lagged number of\nself-determination\nmovements",
  "Past civil war")
modelnames <- c(
  "Base model","SD movements","Past civil war",
  "Random effects","With all checks")

pdf("coefOns.pdf",width=6.1)
OnsetPlot<-multiplot(
  logOnset,logOnset2,logOnset3,logOnset4,logOnset5,intercept=FALSE,
  dodgeHeight=0.5,names=modelnames,
  coefficients=c("logfactions","prevconcessions_l","democracy","kin",
    "yrnocwonset","lagmovements","pastcwonset"))
OnsetPlot + theme_minimal() + ggtitle("Onset models") + 
  scale_color_brewer(palette="Set1",name="") +
  scale_y_discrete(labels=varnames) + 
  theme(axis.text.y = element_text(hjust=0),legend.position="bottom",
    axis.title.y=element_blank(),axis.title.x=element_blank())
dev.off()

# Incidence plot
varnames <- c(
  "Log number of\nself-determination\n factions",
  "Previous concessions","Democracy","Kin in neighboring\nstate",
  "Years since civil\nwar",
  "Lagged number of\nself-determination\nmovements",
  "Past civil war")
modelnames <- c(
  "Base model","SD movements","Past civil war",
  "Random effects","With all checks")

pdf("coefInc.pdf",width=5.75)
IncidencePlot <- multiplot(logIncidence,logIncidence2,logIncidence3,
  logIncidence4,logIncidence5,
  intercept=FALSE,names=modelnames,dodgeHeight=0.5,
  coefficients=c("logfactions","prevconcessions_l","democracy","kin",
  "yrsnocivilwar","lagmovements","pastcw"))
IncidencePlot + theme_minimal() + 
  ggtitle("Incidence models") + 
  scale_color_brewer(palette="Set1",name="") +
  scale_y_discrete(labels=NULL) + 
  theme(axis.text.y=element_text(hjust=0),legend.position="bottom",
    axis.title.y=element_blank(),axis.title.x=element_blank())
dev.off()

###########################
# PREDICTED PROBABILITIES #
###########################

# Sequence of factions
factions.seq <- seq(1,20,0.1)

attach(OnsetClean) # ONSET
# Base model
OnsetX <- cbind(1,log(factions.seq),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrnocwonset),mean(spline1_cwonset),
  mean(spline2_cwonset),mean(spline3_cwonset))
B.Onset <- coef(logOnset)
predprob.Onset <- as.vector(plogis(B.Onset%*%t(OnsetX)))
# lag movements
Onset2X <- cbind(1,log(factions.seq),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrnocwonset),mean(spline1_cwonset),
  mean(spline2_cwonset),mean(spline3_cwonset),mean(lagmovements))
B.Onset2 <- coef(logOnset2)
predprob.Onset2 <- as.vector(plogis(B.Onset2%*%t(Onset2X)))
# Past civil war
Onset3X <- cbind(1,log(factions.seq),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrnocwonset),mean(spline1_cwonset),
  mean(spline2_cwonset),mean(spline3_cwonset),mean(pastcwonset))
B.Onset3 <- coef(logOnset3)
predprob.Onset3 <- as.vector(plogis(B.Onset3%*%t(Onset3X)))
# Random effects
Onset4X <- cbind(1,log(factions.seq),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrnocwonset),mean(spline1_cwonset),
  mean(spline2_cwonset),mean(spline3_cwonset),1)
B.Onset4 <- c(fixef(logOnset4),mean(unlist(ranef(logOnset4,drop=TRUE))))
predprob.Onset4 <- as.vector(plogis(B.Onset4%*%t(Onset4X)))
# All controls
Onset5X <- cbind(1,log(factions.seq),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrnocwonset),mean(spline1_cwonset),
  mean(spline2_cwonset),mean(spline3_cwonset),mean(lagmovements),mean(pastcwonset),1)
B.Onset5 <- c(fixef(logOnset5),mean(unlist(ranef(logOnset5,drop=TRUE))))
predprob.Onset5 <- as.vector(plogis(B.Onset5%*%t(Onset5X)))

detach(OnsetClean)

probs.Onset <- data.frame(cbind(
  factions.seq,predprob.Onset,predprob.Onset2,predprob.Onset3,
  predprob.Onset4,predprob.Onset5))
names(probs.Onset) <- 
  c("factions.seq","base","lagmove","past","re","all")

attach(IncidenceClean) # INCIDENCE
# Base model
IncidenceX <- cbind(1,log(factions.seq),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrsnocivilwar),
  mean(spline1_cw),mean(spline2_cw),mean(spline3_cw))
B.Incidence <- coef(logIncidence)
predprob.Incidence <- as.vector(plogis(B.Incidence%*%t(IncidenceX)))
# lag movements
Incidence2X <- cbind(1,log(factions.seq),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrsnocivilwar),
  mean(spline1_cw),mean(spline2_cw),mean(spline3_cw),mean(lagmovements))
B.Incidence2 <- coef(logIncidence2)
predprob.Incidence2 <- as.vector(plogis(B.Incidence2%*%t(Incidence2X)))
# Past civil war
Incidence3X <- cbind(1,log(factions.seq),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrsnocivilwar),
  mean(spline1_cw),mean(spline2_cw),mean(spline3_cw),mean(pastcw))
B.Incidence3 <- coef(logIncidence3)
predprob.Incidence3 <- as.vector(plogis(B.Incidence3%*%t(Incidence3X)))
# Random effects
Incidence4X <- cbind(1,log(factions.seq),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrsnocivilwar),
  mean(spline1_cw),mean(spline2_cw),mean(spline3_cw),1)
B.Incidence4 <- c(
  fixef(logIncidence4),mean(unlist(ranef(logIncidence4,drop=TRUE))))
predprob.Incidence4 <- as.vector(plogis(B.Incidence4%*%t(Incidence4X)))
# All controls
Incidence5X <- cbind(1,log(factions.seq),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrsnocivilwar),
  mean(spline1_cw),mean(spline2_cw),mean(spline3_cw),
  mean(lagmovements),mean(pastcw),1)
B.Incidence5 <- c(
  fixef(logIncidence5),mean(unlist(ranef(logIncidence5,drop=TRUE))))
predprob.Incidence5 <- as.vector(plogis(B.Incidence5%*%t(Incidence5X)))

detach(IncidenceClean)

probs.Incidence <- data.frame(cbind(
  factions.seq,predprob.Incidence,predprob.Incidence2,
  predprob.Incidence3,predprob.Incidence4,predprob.Incidence5))
names(probs.Incidence) <- 
  c("factions.seq","base","lagmove","past","re","all")


# Plots
probsize <- 0

OnsetPlot <- 
  ggplot(probs.Onset) + 
  geom_line(aes(factions.seq,base,color="Base model"),
    size=probsize) +
  geom_line(aes(factions.seq,lagmove,color="SD movements"),
    size=probsize) +
  geom_line(aes(factions.seq,past,color="Past civil war"),
    size=probsize) +
  geom_line(aes(factions.seq,re,color="Random effects"),
    size=probsize) +
  geom_line(aes(factions.seq,all,color="With all controls"),
    size=probsize)
pdf("probsOnset.pdf",width=5,height=5)
OnsetPlot + theme_minimal() + 
  scale_y_continuous(limits=c(0,0.5)) +
  scale_color_brewer(palette="Set1",name="Model") +
  labs(x="Number of factions",y="") +
  theme(legend.position=c(1,0),legend.justification=c(1,0),
        legend.background = element_rect()) +
  ggtitle("Predicted probability of civil war onset")
dev.off()

  
IncidencePlot <- 
  ggplot(probs.Incidence) + 
  geom_line(aes(factions.seq,base,color="Base model"),
    size=probsize) +
  geom_line(aes(factions.seq,lagmove,color="SD movements"),
    size=probsize) +
  geom_line(aes(factions.seq,past,color="Past civil war"),
    size=probsize) +
  geom_line(aes(factions.seq,re,color="Random effects"),
    size=probsize) +
  geom_line(aes(factions.seq,all,color="With all controls"),
    size=probsize)
pdf("probsIncidence.pdf",width=5,height=5)
IncidencePlot + theme_minimal() + 
  scale_y_continuous(limits=c(0,0.5)) +
  scale_color_brewer(palette="Set1",name="Model") +
  labs(x="Number of factions",y="") +
  theme(legend.position=c(1,0),legend.justification=c(1,0),
        legend.background = element_rect()) +
  ggtitle("Predicted probability of civil war incidence")
dev.off()
