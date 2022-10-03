# 10/02/2022.
# Analysis script for the background motion experiment
# last update: 10/02/22
# The Start ##
# clear the workspace
rm(list=ls()) 
# load libraries in bulk
x<-c("ggpubr", "ggplot2", "multcomp", "pastecs", "tidyr","dplyr", "ggiraph", "ggiraphExtra", "plyr", 
     "covreg", "plot3D", "Hmisc", "corrplot", "psych", "tidyverse", "hrbrthemes", "viridis", "gapminder",
     "ggExtra", "scatterplot3d", "reshape2", "rlang", "plyr", "data.table", "lme4", "magrittr", "fitdistrplus",
     "gridExtra", "statmod", "dotwhisker", "lmerTest")
require(x)
lapply(x, require, character.only = TRUE)
source("~/Documents/GitHub/ger_R_functions/plot_functions.R")
#remove scientific notation in the entire R session
options(scipen = 100)
# locate the data and import it
setwd("/Users/gorkem.er/Desktop/21Projects/background-motion/Data_analysis")
# load bgm data 
bgmdata = read.csv("bgmdata.csv",header=TRUE, quote="\"") 
bgmdata.c = subset(bgmdata, bgmdata$randomTrialsR1C0 == 0)
bgmdata.r = subset(bgmdata, bgmdata$randomTrialsR1C0 == 1)
# load bgmn data 
bgmndata = read.csv("bgmndata.csv",header=TRUE, quote="\"") 
bgmndata.c = subset(bgmndata, bgmndata$randomTrialsR1C0 == 0)
bgmndata.r = subset(bgmndata, bgmndata$randomTrialsR1C0 == 1)
# fgm -> bgm, # note here that "fgmdata" is used to match the naming used in the E1 experiment.
fgmdata <- bgmdata.c 
# FOR RANDOM MOTION ANALYSES bgmdata.r
# fgmdata <- bgmdata.r
number_of_sub <- unique(fgmdata$sub)
# convert response to normed 
fgmdata$responseAR_normed <- bgmndata.c$normedR_indv
# fgmdata$responseAR_normed <- bgmndata.r$normedR_indv # FOR RANDOM ANALYSES
fgmdata$response_error_normed <- bgmndata$normedR_indv - bgmndata$cuedAR
#check the attention for below script!
#fgmdata <- fgmdata[fgmdata$randomTrialsR1C0==1,] ### ATTENTION! HERE TOGGLE ON/OFF random trials! 
# new meta variables
fgmdata$cuedAR <- round(fgmdata$cuedAR, digits = 2)
fgmdata$uncuedAR <- round(fgmdata$uncuedAR, digits = 2)
fgmdata$responseAR <- round(fgmdata$responseAR, digits = 2)
fgmdata$uncuedCat = ifelse(fgmdata$uncuedAR < 0, -1, ifelse(fgmdata$uncuedAR==-0, 0, 1))
fgmdata$uncuedCat <- as.factor(fgmdata$uncuedCat)
fgmdata$cuedCat = ifelse(fgmdata$cuedAR < 0, -1, ifelse(fgmdata$cuedAR==-0, 0, 1))
fgmdata$respAcc <- ifelse( (fgmdata$responseAR > 0 & fgmdata$cuedCat == 1)  | (fgmdata$responseAR < 0 & fgmdata$cuedCat == -1) | (fgmdata$responseAR == 0 & fgmdata$cuedCat == 0), 1, 0)
fgmdata$globalMotion <- as.factor(ifelse(fgmdata$cued_motion_dir == 90 | fgmdata$cued_motion_dir == 270, 1, -1))
# iterate all and get some stats around the data
fgmdata.cleaning <- data.frame(matrix(ncol = 8, nrow = length(number_of_sub)))
colnames(fgmdata.cleaning) <- c("reg_p", "reg_beta", "meanRT", "corr","corr_p", "respAcc","trialN","id")
for (s in 1:length(number_of_sub)){
  tmpdata <- fgmdata[fgmdata$sub == number_of_sub[s],]
  lm_sub <- lm(formula = responseAR ~ cuedAR, data = tmpdata)
  fgmdata.cleaning[s,1] <- round(summary(lm_sub)$coefficients[2,4], digits = 3) #p-value
  fgmdata.cleaning[s,2] <- round(summary(lm_sub)$coefficients[2,1], digits= 3) #estimate
  fgmdata.cleaning[s,3] <- round(mean(tmpdata$rt), digits = 0)/1000
  fgmdata.cleaning[s,4] <- cor.test(tmpdata$cuedAR, tmpdata$responseAR, method = "pearson")$estimate
  fgmdata.cleaning[s,5] <- round(cor.test(tmpdata$cuedAR, tmpdata$responseAR, method = "pearson")$p.value, digits = 2)
  fgmdata.cleaning[s,6] <- mean(tmpdata$respAcc)
  fgmdata.cleaning[s,7] <- nrow(tmpdata)# check the n of each participant
  fgmdata.cleaning[s,8] <- number_of_sub[s]
}
#peak
head(fgmdata.cleaning)
plot(fgmdata.cleaning$reg_beta, fgmdata.cleaning$reg_p)
plot(fgmdata.cleaning$meanRT)
plot(fgmdata.cleaning$corr, fgmdata.cleaning$respAcc) # correlation and accuracy appears to be linked
# possible restrictions
below40P<-fgmdata.cleaning$id[fgmdata.cleaning$respAcc <= 0.40] # 3/74 showed below 0.40 accuracy
below50P<-fgmdata.cleaning$id[fgmdata.cleaning$respAcc <= 0.50] # 3/74 showed below 0.40 accuracy
testToBeCleaned = below50P[12]
fgmdata.cleaning[fgmdata.cleaning$id == testToBeCleaned,]
plot(fgmdata$cuedAR[fgmdata$sub == testToBeCleaned], fgmdata$responseAR[fgmdata$sub == testToBeCleaned])
# those below40P all show bad performance in all metrics (e.g. regression beta, p-value, correlation, reaction time). 
incompletedPeople <- fgmdata.cleaning$id[fgmdata.cleaning$trialN<241]
incompletedPeople
#count globally the trial N
fgmdata.globalTrialN <- data.frame(matrix(ncol = 1, nrow = length(number_of_sub)))
for (s in 1:length(number_of_sub)){
  tmpdata <- fgmdata[fgmdata$sub == number_of_sub[s],]
  fgmdata.globalTrialN[s,1] <- nrow(tmpdata)# check the n of each participant
  fgmdata.globalTrialN[s,2] <- number_of_sub[s]
}
colnames(fgmdata.globalTrialN) <- c("trialN","id")
head(fgmdata.globalTrialN)
incompletedPeople_global <- fgmdata.globalTrialN$id[fgmdata.globalTrialN$trialN<241]
plot(fgmdata.globalTrialN[,2], fgmdata.globalTrialN[,1])
abline(h = 241)
# check coherent and global trial N people
incompletedPeople
incompletedPeople_global
# for random analyses, I need to remove those people using the coherent analysis
# incompletedPeople <- c(19117, 19142, 19242) # WARNING THIS IS FOR RANDOM ANALYSIS ONLY
#### CLEANED DATA ####
# removing incomplete experiment people
fgmdata <- fgmdata[!( (fgmdata$sub %in% incompletedPeople)),]
# check, also remove id 18667 if it is there
18667 %in% unique(fgmdata$sub) 
fgmdata <- fgmdata[!(fgmdata$sub == 18667),] # because 0.53 corr
below50P[1:length(below50P)] %in% unique(fgmdata$sub) # many of the below50 people are included in the analysis, just fyi
# some descriptive plots #
findSameAR_trials <- summarySE(fgmdata, measurevar="responseAR", groupvars=c("cuedAR", "identicalShapesI1D0", "sub"))
sameAR_trials <- findSameAR_trials[findSameAR_trials$identicalShapesI1D0==1,]
sameAR_trials #so many less-same trials per individuals
par(mfrow = c(2,2))
plot(findSameAR_trials$cuedAR[findSameAR_trials$identicalShapesI1D0==1], findSameAR_trials$responseAR[findSameAR_trials$identicalShapesI1D0==1], xlab = "Cued AR", ylab = "Response AR", main = "Same-AR Trials")
abline(coef = c(0,1),col="red", lwd=3, lty=2)
plot(findSameAR_trials$cuedAR[findSameAR_trials$identicalShapesI1D0==0], findSameAR_trials$responseAR[findSameAR_trials$identicalShapesI1D0==0], xlab = "Cued AR", ylab = "Response AR",main = "Different-AR Trials", ylim = c(-0.6, 0.6))
abline(coef = c(0,1),col="red", lwd=3, lty=2)
#uncuedAR
findSameAR_trials_uncued <- summarySE(fgmdata, measurevar="responseAR", groupvars=c("uncuedAR", "identicalShapesI1D0", "sub"))
plot(findSameAR_trials_uncued$uncuedAR[findSameAR_trials_uncued$identicalShapesI1D0==1], findSameAR_trials_uncued$responseAR[findSameAR_trials_uncued$identicalShapesI1D0==1], xlab = "Uncued AR", ylab = "Response AR", main = "Same-AR Trials")
abline(coef = c(0,1),col="red", lwd=3, lty=2)
plot(findSameAR_trials_uncued$uncuedAR[findSameAR_trials_uncued$identicalShapesI1D0==0], findSameAR_trials_uncued$responseAR[findSameAR_trials_uncued$identicalShapesI1D0==0], xlab = "Uncued AR", ylab = "Response AR", main = "Different-AR Trials", ylim = c(-0.6, 0.6))
abline(coef = c(0,1),col="red", lwd=3, lty=2)
abline(coef = c(0,0),col="red", lwd=3, lty=2)
par(mfrow = c(1,1))
#### FINAL ADDING TO THIS SCRIPT #### 
# main analyses and writing out the regression formulas #
library(equatiomatic)
# 1) regression plot stats
fullModel <- lmer(response_error ~ uncuedAR * sameDirection1S0D + (1 | sub), data = fgmdata, REML = FALSE)
summary(fullModel)
extract_eq(fullModel, wrap = TRUE, terms_per_line = 2)
# 2) uncued beta coeff stats
number_of_sub <- unique(fgmdata$sub)
tmpdata <- aggregate(response_error~ uncuedAR + sub + sameDirection1S0D, fgmdata, mean)
fgmdata.indv_beta <- data.frame(matrix(ncol = 3, nrow = length(number_of_sub)))
for (r in 1:length(number_of_sub)){ 
  tmpdata_sub <- tmpdata[tmpdata$sub==number_of_sub[r],]
  #run a regression model on individual sub
  lm_sub_diff <- lm(response_error ~ uncuedAR, data = tmpdata_sub[tmpdata_sub$sameDirection1S0D==0,])
  lm_beta_diff <- summary(lm_sub_diff)$coefficients[2]
  lm_sub_same <- lm(response_error ~ uncuedAR, data = tmpdata_sub[tmpdata_sub$sameDirection1S0D==1,])
  lm_beta_same <- summary(lm_sub_same)$coefficients[2]
  fgmdata.indv_beta[r,1] <- lm_beta_diff
  fgmdata.indv_beta[r,2] <- lm_beta_same
  fgmdata.indv_beta[r,3] = number_of_sub[r]
}
head(fgmdata.indv_beta)
plot(fgmdata.indv_beta$X1, fgmdata.indv_beta$X2)
abline(c(0,1)) # more points lie left of the abline, same has higher response errors
#long to wide format
meltData <- melt(fgmdata.indv_beta[1:2])
my_comparisons = list(c("X1","X2"))
gglinePlot <- ggline(meltData, x = "variable", y = "value", 
                     add = c("mean_ci", "jitter"), palette = "jco")+ 
  stat_compare_means(paired = TRUE, comparisons = my_comparisons)+
  stat_compare_means(label.y = 0.3)
gglinePlot
summary(gglinePlot)
compare_means(value ~ variable, data = meltData, paired = TRUE,  method = "t.test")# alternative = "greater", method = "t.test"
t.test(meltData$value[meltData$variable == "X1"], meltData$value[meltData$variable == "X2"], paired = T)
# 3) global organization stats
fgmdata$global_org <- fgmdata$global_org1W0B # changing naming for compatability with the E1 results
globalOrg_model <- lmer(response_error ~ uncuedAR * sameDirection1S0D * global_org  + (1 | sub), data = fgmdata,  REML = FALSE)
summary(globalOrg_model)
extract_eq(globalOrg_model, wrap = TRUE, terms_per_line = 2)
# doing it same with the ### NORMED RESPONSE AR !! ####
# 1) regression
normedResponseAR_model <- lmer(responseAR_normed ~ uncuedAR * sameDirection1S0D_R2 + (1 | sub), data = fgmdata)
summary(normedResponseAR_model)
extract_eq(fullModel, wrap = TRUE, terms_per_line = 2)
# 2) global org
globalOrg_model_responseAR <- lmer(responseAR_normed ~ uncuedAR * sameDirection1S0D * global_org  + (1 | sub), data = fgmdata,  REML = FALSE)
summary(globalOrg_model_responseAR)
extract_eq(globalOrg_model_responseAR, wrap = TRUE, terms_per_line = 2)
#summary(lmer(responseAR_normed ~ uncuedAR + uncuedAR:sameDirection1S0D_R2 * global_org1W0B + (1 | sub) + (1 | sub:uncuedAR) + (1 | sub:global_org1W0B), data = fgmdata, REML = FALSE))
# 3) beta plot response normed
number_of_sub <- unique(fgmdata$sub)
tmpdata <- aggregate(responseAR_normed~ uncuedAR + sub + sameDirection1S0D, fgmdata, mean)
fgmdata.indv_beta_response_normed <- data.frame(matrix(ncol = 3, nrow = length(number_of_sub)))
for (r in 1:length(number_of_sub)){ 
  tmpdata_sub <- tmpdata[tmpdata$sub==number_of_sub[r],]
  #run a regression model on individual sub
  lm_sub_diff <- lm(responseAR_normed ~ uncuedAR, data = tmpdata_sub[tmpdata_sub$sameDirection1S0D==0,])
  lm_beta_diff <- summary(lm_sub_diff)$coefficients[2]
  lm_sub_same <- lm(responseAR_normed ~ uncuedAR, data = tmpdata_sub[tmpdata_sub$sameDirection1S0D==1,])
  lm_beta_same <- summary(lm_sub_same)$coefficients[2]
  fgmdata.indv_beta_response_normed[r,1] <- lm_beta_diff
  fgmdata.indv_beta_response_normed[r,2] <- lm_beta_same
  fgmdata.indv_beta_response_normed[r,3] = number_of_sub[r]
}
head(fgmdata.indv_beta_response_normed)
plot(fgmdata.indv_beta_response_normed$X1, fgmdata.indv_beta_response_normed$X2)
abline(c(0,1)) # more points lie left of the abline, same has higher response errors
#### MELT DATA ####
#long to wide format
meltData <- melt(fgmdata.indv_beta_response_normed[1:2])
head(meltData)
my_comparisons = list(c("X1","X2"))
beta_plot_response_normed <- ggline(meltData, x = "variable", y = "value",
                                    add = c("mean_ci", "jitter"), add.params = list(color ="variable", size = 2, alpha = 0.5))+
  stat_compare_means(paired= TRUE, comparisons = my_comparisons)
beta_plot_response_normed
compare_means(value ~ variable, data = meltData, paired = TRUE,  method = "t.test")# alternative = "greater", method = "t.test"
t.test(meltData$value[meltData$variable == "X1"], meltData$value[meltData$variable == "X2"], paired = T)
# random motion
# toggle off things at the top
# fgm -> bgm
fgmdata <- bgmdata.r
fgmdata$responseAR_normed <- bgmndata.r$normedR_indv
fgmdata <- fgmdata[fgmdata$randomTrialsR1C0==1,] #random trials! 
# iterate all and get the beta coefficient
fgmdata.cleaning <- data.frame(matrix(ncol = 8, nrow = length(number_of_sub)))
colnames(fgmdata.cleaning) <- c("reg_p", "reg_beta", "meanRT", "corr","corr_p", "respAcc","trialN","id")
for (s in 1:length(number_of_sub)){
  tmpdata <- fgmdata[fgmdata$sub == number_of_sub[s],]
  lm_sub <- lm(formula = responseAR ~ cuedAR, data = tmpdata)
  fgmdata.cleaning[s,1] <- round(summary(lm_sub)$coefficients[2,4], digits = 3) #p-value
  fgmdata.cleaning[s,2] <- round(summary(lm_sub)$coefficients[2,1], digits= 3) #estimate
  fgmdata.cleaning[s,3] <- round(mean(tmpdata$rt), digits = 0)/1000
  fgmdata.cleaning[s,4] <- cor.test(tmpdata$cuedAR, tmpdata$responseAR, method = "pearson")$estimate
  fgmdata.cleaning[s,5] <- round(cor.test(tmpdata$cuedAR, tmpdata$responseAR, method = "pearson")$p.value, digits = 2)
  fgmdata.cleaning[s,6] <- mean(tmpdata$respAcc)
  fgmdata.cleaning[s,7] <- nrow(tmpdata)# check the n of each participant
  fgmdata.cleaning[s,8] <- number_of_sub[s]
}
head(fgmdata.cleaning)
# 21 agu, I found out that 19117 19142 19242 are incomplete, also remove id 18667
incompletedPeople <- c(19117, 19142, 19242, 18667) # WARNING THIS IS FOR RANDOM ANALYSIS ONLY
incompletedPeople
fgmdata <- fgmdata[!( (fgmdata$sub %in% incompletedPeople)),]
# random analysis start
fgmdata$sameDirection1S0D_R2 <- as.factor(fgmdata$sameDirection1S0D_R2)
# response error x uncued AR 
summary(lmer(response_error ~ uncuedAR + (1 | sub) + (1 | sub:uncuedAR), data = fgmdata, REML = FALSE))
# global org
summary(lmer(response_error ~ uncuedAR * global_org1W0B + (1 | sub) + (1 | sub:uncuedAR) + (1 | sub:global_org1W0B), data = fgmdata, REML = FALSE))
# normed response
summary(lmer(responseAR_normed ~ uncuedAR + (1 | sub) + (1 | sub:uncuedAR), data = fgmdata, REML = FALSE))
# global org
summary(lmer(responseAR_normed ~ uncuedAR * global_org1W0B + (1 | sub) + (1 | sub:uncuedAR) + (1 | sub:global_org1W0B), data = fgmdata, REML = FALSE))
# THE END. # 

