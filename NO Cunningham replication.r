# NOTES: Look at residual deviance instead of just residuals
# Can still look at temporal or spatial clustering even with a logit

library(foreign)
library(sandwich)
library(ggplot2)
library(lmtest)
library(boot)
library(MASS)
library(stargazer)
library(RColorBrewer)

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
  "spline1_cwonset","spline2_cwonset","spline3_cwonset","spline1_cw",
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
  kgcid,country,group,year,acdcivilwar1,logfactions,prevconcessions_l,
  democracy,kin,yrsnocivilwar,spline1_cw,spline2_cw,spline3_cw))
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

# Check models
summary(logOnsetTest)
summary(logIncidenceTest)

varnames <- c(
  "Log SD factions","Previous concessions","Democracy","Kin",
  "Years since war","Years since onset")
se <- list(logOnsetTest[,"Std. Error"],logIncidenceTest[,"Std. Error"])
p <- list(logOnsetTest[,"Pr(>|z|)"],logIncidenceTest[,"Pr(>|z|)"])
		
# REPLICATE TABLE
sink("RepTable.tex")		
stargazer(logOnset,logIncidence,se=se,p=p,omit=c(6:8,10:12),
  covariate.labels=varnames,model.numbers=FALSE,
  dep.var.labels=c("Civil war onset","Civil war incidence"),
  digits=2,align=TRUE,style="ajps",float=FALSE)
sink()

# REPLICATE PREDICTED PROBABILITIES

# Sequence of factions
factions.seq <- seq(1,20,0.1)

# Set up design matrix (yes previous concessions, no democracy, yes kin)
attach(OnsetClean)
OnsetX <- cbind(1,log(factions.seq),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrnocwonset),mean(spline1_cwonset),
  mean(spline2_cwonset),mean(spline3_cwonset))
detach(OnsetClean)

attach(IncidenceClean)
IncidenceX <- cbind(1,log(factions.seq),median(prevconcessions_l),
  median(democracy),median(kin),mean(yrsnocivilwar),mean(spline1_cw),
  mean(spline2_cw),mean(spline3_cw))
detach(IncidenceClean)

B.Onset <- mvrnorm(10000, coef(logOnset), logOnsetClust)
predprob.Onset <- plogis(B.Onset%*%t(OnsetX))
predprob.Onset <- t(apply(predprob.Onset,2,quantile,c(0.025,0.925)))
predprob.Onset <- data.frame(cbind(
  factions.seq,predprob.Onset,apply(predprob.Onset,1,mean)))
colnames(predprob.Onset) <- c("factions.seq","ci.low","ci.high","mean")

B.Incidence <- mvrnorm(10000, coef(logIncidence), logIncidenceClust)
predprob.Incidence <- plogis(B.Incidence%*%t(IncidenceX))
predprob.Incidence <- 
  t(apply(predprob.Incidence,2,quantile,c(0.025,0.925)))
predprob.Incidence <- data.frame(cbind(
  factions.seq,predprob.Incidence,apply(predprob.Incidence,1,mean)))
colnames(predprob.Incidence) <- c(
  "factions.seq","ci.low","ci.high","mean")

# Plot simulated confidence intervals
OnsetPlot <- 
  ggplot(predprob.Onset, aes(factions.seq,mean)) + geom_line() + 
  geom_ribbon(aes(ymin=ci.low,ymax=ci.high),alpha=0.3,fill="blue")
IncidencePlot <- 
  ggplot(predprob.Incidence, aes(factions.seq,mean)) + geom_line() + 
  geom_ribbon(aes(ymin=ci.low,ymax=ci.high),alpha=0.3,fill="blue")

pdf("predprobOnset.pdf",width=5,height=5)
OnsetPlot + theme_minimal() + 
  scale_x_continuous(limits=c(0,max(factions.seq)),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,0.6)) +
  labs(x="Number of factions",y="") +
  ggtitle("Predicted probability of civil war onset")
dev.off()
pdf("predprobIncidence.pdf",width=5,height=5)
IncidencePlot + theme_minimal() + 
  scale_x_continuous(limits=c(0,max(factions.seq)),expand=c(0,0)) +
  scale_y_continuous(limits=c(0,0.6)) +
  labs(x="Number of factions",y="") +
  ggtitle("Predicted probability of civil war incidence")
dev.off()