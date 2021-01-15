library(readxl)
library(ggplot2)
library(lmerTest)
library(dplyr)
library(tidyr)
library(PMCMRplus)
library(repolr)
library(lm.beta)
library(ggbeeswarm)
library(tidyverse)
library(anchors)


## SETUP ##

setwd("/Users/nikol/Documents/LEQ PRoject/Paper 2-2020/May2020/ReviewerComments_11-18/Cog_Stuff/")
#LEQ<-read_xlsx("./LEQ_Data_10.1.19(3).xlsx", sheet="Sheet1")
LEQ<-read_xlsx("./LEQ_Data_12-9-2020.xlsx", sheet="INQuery")
LEQ$INDDID<-as.factor(LEQ$INDDID)
LEQ$Group<-as.factor(LEQ$Group)
LEQ$Sex<-as.factor(LEQ$Sex)
LEQ$Visit<-as.factor(LEQ$Visit)
LEQ$Long<-as.factor(LEQ$Long)

##Normality Checking##
hist(LEQ$Descriptor_Total_T)
shapiro.test(LEQ$Descriptor_Total_T)

hist(LEQ$BehavioralScale)
shapiro.test(LEQ$BehavioralScale)

hist(LEQ$BensonCopy)
shapiro.test(LEQ$BensonCopy)

hist(LEQ$BensonRecall)
shapiro.test(LEQ$BensonRecall)

hist(LEQ$FWordsCorrect)
shapiro.test(LEQ$FWordsCorrect)

hist(LEQ$BNT)
shapiro.test(LEQ$BNT)

hist(LEQ$ExecutiveScale)
shapiro.test(LEQ$ExecutiveScale) 

hist(LEQ$LanguageScale)
shapiro.test(LEQ$LanguageScale) 

hist(LEQ$VisualScale)
shapiro.test(LEQ$VisualScale) 

hist(LEQ$TotalMemoryScale)
shapiro.test(LEQ$TotalMemoryScale) 

##Tukey's ladder -- define function##
DataTransformation = function(var){
  df<-as.data.frame(var) %>% transmute("Rank" = rank(var, ties.method = "average", na.last = "keep"), "lambda=3" = var^3, "lambda=2" = var^2, "lambda=1" = var, "lambda=0.5" = var^(1/2), "lambda=0" = log(var), "lambda=-0.5" = -var^(-1/2), "lambda=-1" = -var^(-1), "lambda=-2" = -var^(-2), "lambda=-3" = -var^(-3))
  df<-anchors::replace.value(data=df, names = c("lambda=3", "lambda=2", "lambda=1", "lambda=0.5", "lambda=0", "lambda=-0.5", "lambda=-1", "lambda=-2", "lambda=-3"), from = "-Inf", to = as.double(NA))
  print(lapply(df, shapiro.test))
  print(gather(df, key = "Power", value = "Value", -Rank) %>% mutate(Power = factor(Power, levels = c("lambda=3", "lambda=2", "lambda=1", "lambda=0.5", "lambda=0", "lambda=-0.5", "lambda=-1", "lambda=-2", "lambda=-3"))) %>% ggplot(aes(x=0, y=Value)) + geom_boxplot(outlier.shape = NA) + geom_beeswarm(shape=21, aes(fill = Rank)) + facet_wrap(~Power, scales = "free_y") + labs(x="", y="", title = "Tukey's ladder of powers y=x^lambda"))
}

##Rank Transform then Tukey's Ladder Raw Variables##
DataTransformation(LEQ$Descriptor_Total_T) 
DataTransformation(LEQ$BehavioralScale)
DataTransformation(LEQ$BNT) 
DataTransformation(LEQ$FWordsCorrect) 
DataTransformation(LEQ$BensonCopy) 
DataTransformation(LEQ$BensonRecall) 
DataTransformation(LEQ$ExecutiveScale) 
DataTransformation(LEQ$LanguageScale) 
DataTransformation(LEQ$VisualScale) 
DataTransformation(LEQ$TotalMemoryScale) 
DataTransformation(LEQ$PBACTotal) 


##Format CrossSectional + Longitudinal Cohorts##
LEQ_t1<-subset(LEQ, Visit=="1")
LEQ_Long<-subset(LEQ, Long=="Y")
LEQ_Long_t1<-subset(LEQ_Long, Visit=="1")


#Summary
summary(LEQ$Group)
summary(LEQ_t1$AgeatTest_CDR)
summary(LEQ$Sex)
summary(LEQ$Education)
sd(LEQ$Education,na.rm=TRUE)
summary(LEQ$Overall_DD_MRI)
sd(LEQ$Overall_DD_MRI,na.rm=TRUE)
summary(LEQ$FTLDCDRSum)
sd(LEQ$FTLDCDRSum,na.rm=TRUE)



# Original Individual Variable Models -------------------------------------

## Social Behavior Observer (SBO) Checklist - Descriptor Total ##

#Cross Sectional Late LEQ
summary(lm(LEQ_t1$Descriptor_Total_T~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_SBO + LEQ_t1$DD_SBO_YR, na.action = "na.exclude"))
ggplot(data=LEQ_t1, aes(y=LEQ_t1$Descriptor_Total_T, x=LEQ_t1$Late_LEQ, color=Group)) + geom_line(aes(group=INDDID), linetype=3, size=1) + geom_point(aes(size=1)) + guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2) + theme_classic(base_size = 17)
lm.beta::lm.beta(lm(LEQ_t1$Descriptor_Total_T~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_SBO + LEQ_t1$DD_SBO_YR, na.action = "na.exclude"))
plot(lm(LEQ_t1$Descriptor_Total_T~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_SBO + LEQ_t1$DD_SBO_YR, na.action = "na.exclude"))
shapiro.test((lm(LEQ_t1$Descriptor_Total_T~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_SBO + LEQ_t1$DD_SBO_YR, na.action = "na.exclude")$residuals))
SBOCS <- summary(lm(LEQ_t1$Descriptor_Total_T~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_SBO + LEQ_t1$DD_SBO_YR, na.action = "na.exclude"))$r.squared
SBOCS / (1-SBOCS)

#Longitudinal
ggplot(data=LEQ_Long, aes(y=Descriptor_Total_T, x=DD_SBO_YR)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)

SBO_DTslope<-lmer(Descriptor_Total_T~(DD_SBO_YR|INDDID)+ SexFormatted + AgeatTestSBO_base, data=LEQ_Long, na.action = "na.exclude", control=lmerControl(check.nobs.vs.nRE="ignore")) 
summary(SBO_DTslope)
SBO_DTslopes<-coef(SBO_DTslope)$INDDID
SBO_DTslopes$DD_MRI
SBO_DTslopes<-setNames(cbind(rownames(SBO_DTslopes),SBO_DTslopes,row.names=NULL),c("INDDID", "SBO_DTSlopes","SBO_DTIntercept","NA1","NA2"))
SBO_DTkeep<-SBO_DTslopes
LEQ_Long_slopes<-merge(LEQ_Long_t1, SBO_DTkeep, by="INDDID",all.x = TRUE)

##Late##
summary(lm(LEQ_Long_slopes$SBO_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude"))
ggplot(data=LEQ_Long_slopes, aes(y=SBO_DTSlopes, x=LEQ_Long_t1$Late_LEQ)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20) + labs(x="Late-Life LEQ", y="Slope of SBOC Descriptor Score")
lm.beta::lm.beta(lm(LEQ_Long_slopes$SBO_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude"))
plot(lm(LEQ_Long_slopes$SBO_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude"))
shapiro.test(lm(LEQ_Long_slopes$SBO_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude")$residuals)
SBOLong <- summary(lm(LEQ_Long_slopes$SBO_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude"))$r.squared
SBOLong / (1-SBOLong)


## PBAC Behavior ##

#Cross Sectional Late LEQ
summary(lm(LEQ_t1$BehavioralScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$DD_PBAC_YR))
ggplot(data=LEQ_t1, aes(y=LEQ_t1$BehavioralScale, x=LEQ_t1$Late_LEQ)) + geom_line(aes(group=INDDID), linetype=3, size=1) + geom_point(aes(size=1)) + guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2) + theme_classic(base_size = 17) + labs(x="Late-Life LEQ", y="PBAC Behavior Score")
lm.beta::lm.beta(lm(LEQ_t1$BehavioralScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$DD_PBAC_YR))
plot(lm(LEQ_t1$BehavioralScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$DD_PBAC_YR))
shapiro.test(lm(LEQ_t1$BehavioralScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$DD_PBAC_YR)$residuals)
PBACbeCS <- summary(lm(LEQ_t1$BehavioralScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$DD_PBAC_YR))$r.squared
PBACbeCS / (1-PBACbeCS)


#Longitudinal
ggplot(data=LEQ_Long, aes(y=BehavioralScale, x=DD_PBAC_YR)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)

PBACBehslope<-lmer(BehavioralScale~(DD_PBAC_YR|INDDID)+ AgeatTestPBAC_base + SexFormatted, data=LEQ_Long, control=lmerControl(check.nobs.vs.nRE="ignore")) 
summary(PBACBehslope)
PBACBehslopes<-coef(PBACBehslope)$INDDID
PBACBehslopes$DD_MRI
PBACBehslopes<-setNames(cbind(rownames(PBACBehslopes),PBACBehslopes,row.names=NULL),c("INDDID", "PBACBehslopes","PBACBehIntercept","NA1","NA2"))
PBACBehkeep<-PBACBehslopes
LEQ_Long_slopes<-merge(LEQ_Long_t1, PBACBehkeep, by="INDDID",all.x = TRUE)

##Late##
summary(lm(LEQ_Long_slopes$PBACBehslopes~LEQ_Long_t1$Late_LEQ))
ggplot(data=LEQ_Long_slopes, aes(y=PBACBehslopes, x=LEQ_Long_t1$Late_LEQ)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)
lm.beta::lm.beta(lm(LEQ_Long_slopes$PBACBehslopes~LEQ_Long_t1$Late_LEQ))
plot(lm(LEQ_Long_slopes$PBACBehslopes~LEQ_Long_t1$Late_LEQ))
shapiro.test(lm(LEQ_Long_slopes$PBACBehslopes~LEQ_Long_t1$Late_LEQ)$residuals)
PBACbeLongR2 <- summary(lm(LEQ_Long_slopes$PBACBehslopes~LEQ_Long_t1$Late_LEQ))$r.squared
PBACbeLongR2 / (1-PBACbeLongR2)


## F Words ##

#Cross Sectional Late LEQ
summary(lm(LEQ_t1$FWordsCorrect~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_VF + LEQ_t1$DD_Fletter_YR))
lm.beta::lm.beta(lm(LEQ_t1$FWordsCorrect~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_VF + LEQ_t1$DD_Fletter_YR))
plot(lm(LEQ_t1$FWordsCorrect~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_VF + LEQ_t1$DD_Fletter_YR))
shapiro.test(lm(LEQ_t1$FWordsCorrect~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_VF + LEQ_t1$DD_Fletter_YR)$residuals)




#Longitudinal
ggplot(data=LEQ_Long, aes(y=FWordsCorrect, x=DD_Fletter_YR)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)+labs(x="Disease Duration (years)", y="F Words")

Fslope<-lmer(FWordsCorrect~(DD_Fletter_YR|INDDID)+ AgeatTestFWords_base + SexFormatted, data=LEQ_Long, na.action = "na.exclude", control=lmerControl(check.nobs.vs.nRE="ignore"))
summary(Fslope)
Fslopes<-coef(Fslope)$INDDID
Fslopes$DD_MRI
Fslopes<-setNames(cbind(rownames(Fslopes),Fslopes,row.names=NULL),c("INDDID", "Fslopes","FIntercept","NA1","NA2"))
Fkeep<-Fslopes
LEQ_Long_slopes<-merge(LEQ_Long_t1, Fkeep, by="INDDID",all.x = TRUE)

#Late#
summary(lm(LEQ_Long_slopes$Fslopes~LEQ_Long_t1$Late_LEQ))
lm.beta::lm.beta(lm(LEQ_Long_slopes$Fslopes~LEQ_Long_t1$Late_LEQ))
plot(lm(LEQ_Long_slopes$Fslopes~LEQ_Long_t1$Late_LEQ))
shapiro.test(lm(LEQ_Long_slopes$Fslopes~LEQ_Long_t1$Late_LEQ)$residuals)



## BNT ##

#Cross Sectional Late LEQ
summary(lm(LEQ_t1$BNT~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_BNT + LEQ_t1$DD_BNT_YR))
lm.beta::lm.beta(lm(LEQ_t1$BNT~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_BNT + LEQ_t1$DD_BNT_YR))
plot(lm(LEQ_t1$BNT~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_BNT + LEQ_t1$DD_BNT_YR))
shapiro.test(lm(LEQ_t1$BNT~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_BNT + LEQ_t1$DD_BNT_YR)$residuals)



#Longitudinal
ggplot(data=LEQ_Long, aes(y=BNT, x=DD_BNT_YR)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)+labs(x="Disease Duration (years)", y="BNT")

BNTslope<-lmer(BNT~(DD_BNT_YR|INDDID)+ AgeatTestBNT_base + SexFormatted, data=LEQ_Long, na.action = "na.exclude", control=lmerControl(check.nobs.vs.nRE="ignore")) 
summary(BNTslope)
BNTslopes<-coef(BNTslope)$INDDID
BNTslopes$DD_MRI
BNTslopes<-setNames(cbind(rownames(BNTslopes),BNTslopes,row.names=NULL),c("INDDID", "BNTslopes","BNTIntercept","NA1","NA2"))
BNTkeep<-BNTslopes
LEQ_Long_slopes<-merge(LEQ_Long_t1, BNTkeep, by="INDDID",all.x = TRUE)

##Late##
summary(lm(LEQ_Long_slopes$BNTslopes~LEQ_Long_t1$Late_LEQ))
lm.beta::lm.beta(lm(LEQ_Long_slopes$BNTslopes~LEQ_Long_t1$Late_LEQ))
plot(lm(LEQ_Long_slopes$BNTslopes~LEQ_Long_t1$Late_LEQ))
shapiro.test(lm(LEQ_Long_slopes$BNTslopes~LEQ_Long_t1$Late_LEQ)$residuals)



##Benson##

##Cross Sectional Copy+Recall Late LEQ##
summary(lm(LEQ_t1$BensonCopy~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_Benson + LEQ_t1$DD_BensonC_YR, na.action="na.exclude"))
lm.beta::lm.beta(lm(LEQ_t1$BensonCopy~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_Benson + LEQ_t1$DD_BensonC_YR, na.action="na.exclude"))
plot(lm(LEQ_t1$BensonCopy~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_Benson + LEQ_t1$DD_BensonC_YR, na.action="na.exclude"))
shapiro.test(lm(LEQ_t1$BensonCopy~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_Benson + LEQ_t1$DD_BensonC_YR, na.action="na.exclude")$residuals)



summary(lm(LEQ_t1$BensonRecall~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_Benson + LEQ_t1$DD_BensonR_YR, na.action="na.exclude"))
lm.beta::lm.beta(lm(LEQ_t1$BensonRecall~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_Benson + LEQ_t1$DD_BensonR_YR, na.action="na.exclude"))
plot(lm(LEQ_t1$BensonRecall~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_Benson + LEQ_t1$DD_BensonR_YR, na.action="na.exclude"))
shapiro.test(lm(LEQ_t1$BensonRecall~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_Benson + LEQ_t1$DD_BensonR_YR, na.action="na.exclude")$residuals)



#Longitudinal
ggplot(data=LEQ_Long, aes(y=BensonCopy, x=DD_BensonC_YR)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)+labs(x="Disease Duration (years)", y="Benson Copy")
ggplot(data=LEQ_Long, aes(y=BensonRecall, x=DD_BensonR_YR)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)+labs(x="Disease Duration (years)", y="Benson Recall")

##Copy##
BensonCopyslope<-lmer(BensonCopy~(DD_BensonC_YR|INDDID)+ AgeatTestBCopy_base + SexFormatted, data=LEQ_Long, na.action = "na.exclude", control=lmerControl(check.nobs.vs.nRE="ignore"), REML = FALSE) 
summary(BensonCopyslope)
BensonCopyslopes<-coef(BensonCopyslope)$INDDID
BensonCopyslopes$DD_MRI
BensonCopyslopes<-setNames(cbind(rownames(BensonCopyslopes),BensonCopyslopes,row.names=NULL),c("INDDID", "BensonCopySlopes","BensonCopyIntercept","NA1","NA2"))
BensonCopykeep<-BensonCopyslopes
LEQ_Long_slopes<-merge(LEQ_Long_t1, BensonCopykeep, by="INDDID",all.x = TRUE)

##Late##
summary(lm(LEQ_Long_slopes$BensonCopySlopes~LEQ_Long_t1$Late_LEQ))
lm.beta::lm.beta(lm(LEQ_Long_slopes$BensonCopySlopes~LEQ_Long_t1$Late_LEQ))
plot(lm(LEQ_Long_slopes$BensonCopySlopes~LEQ_Long_t1$Late_LEQ))
shapiro.test(lm(LEQ_Long_slopes$BensonCopySlopes~LEQ_Long_t1$Late_LEQ)$residuals)




##Recall##
BensonRecallslope<-lmer(BensonRecall~(DD_BensonR_YR|INDDID)+ AgeatTestBRecall_base + SexFormatted, data=LEQ_Long, na.action = "na.exclude", control=lmerControl(check.nobs.vs.nRE="ignore")) 
summary(BensonRecallslope)
BensonRecallslopes<-coef(BensonRecallslope)$INDDID
BensonRecallslopes$DD_MRI
BensonRecallslopes<-setNames(cbind(rownames(BensonRecallslopes),BensonRecallslopes,row.names=NULL),c("INDDID", "BensonRecallSlopes","BensonRecallIntercept","NA1","NA2"))
BensonRecallkeep<-BensonRecallslopes
LEQ_Long_slopes<-merge(LEQ_Long_t1, BensonRecallkeep, by="INDDID",all.x = TRUE)

#Late#
summary(lm(LEQ_Long_slopes$BensonRecallSlopes~LEQ_Long_t1$Late_LEQ))
lm.beta::lm.beta(lm(LEQ_Long_slopes$BensonRecallSlopes~LEQ_Long_t1$Late_LEQ))
plot(lm(LEQ_Long_slopes$BensonRecallSlopes~LEQ_Long_t1$Late_LEQ))
shapiro.test(lm(LEQ_Long_slopes$BensonRecallSlopes~LEQ_Long_t1$Late_LEQ)$residuals)









# PBAC Global Scales ------------------------------------------------------

##PBAC Executive## 

#Cross Sectional Late LEQ
summary(lm(LEQ_t1$ExecutiveScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$DD_PBAC_YR, na.action = "na.exclude"))
ggplot(data=LEQ_t1, aes(y=LEQ_t1$ExecutiveScale, x=LEQ_t1$Late_LEQ, color=Group)) + geom_line(aes(group=INDDID), linetype=3, size=1) + geom_point(aes(size=1)) + guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2) + theme_classic(base_size = 17)
lm.beta::lm.beta(lm(LEQ_t1$ExecutiveScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$DD_PBAC_YR, na.action = "na.exclude"))
plot(lm(LEQ_t1$ExecutiveScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$DD_PBAC_YR, na.action = "na.exclude"))
shapiro.test((lm(LEQ_t1$ExecutiveScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$DD_PBAC_YR, na.action = "na.exclude")$residuals))
PBACExecCS <- summary(lm(LEQ_t1$ExecutiveScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$DD_PBAC_YR, na.action = "na.exclude"))$r.squared
PBACExecCS / (1-PBACExecCS)

#Longitudinal
ggplot(data=LEQ_Long, aes(y=ExecutiveScale, x=DD_PBAC_YR)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)

PBAC_DTslope<-lmer(ExecutiveScale~(DD_PBAC_YR|INDDID)+ SexFormatted + AgeatTestPBAC_base, data=LEQ_Long, na.action = "na.exclude", control=lmerControl(check.nobs.vs.nRE="ignore")) 
summary(PBAC_DTslope)
PBAC_DTslopes<-coef(PBAC_DTslope)$INDDID
PBAC_DTslopes$DD_MRI
PBAC_DTslopes<-setNames(cbind(rownames(PBAC_DTslopes),PBAC_DTslopes,row.names=NULL),c("INDDID", "PBAC_DTSlopes","PBAC_DTIntercept","NA1","NA2"))
PBAC_DTkeep<-PBAC_DTslopes
LEQ_Long_slopes<-merge(LEQ_Long_t1, PBAC_DTkeep, by="INDDID",all.x = TRUE)

##Late##
summary(lm(LEQ_Long_slopes$PBAC_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude"))
ggplot(data=LEQ_Long_slopes, aes(y=PBAC_DTSlopes, x=LEQ_Long_t1$Late_LEQ)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20) + labs(x="Late-Life LEQ", y="Slope of PBACC Descriptor Score")
lm.beta::lm.beta(lm(LEQ_Long_slopes$PBAC_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude"))
plot(lm(LEQ_Long_slopes$PBAC_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude"))
shapiro.test(lm(LEQ_Long_slopes$PBAC_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude")$residuals)
PBACExecLong <- summary(lm(LEQ_Long_slopes$PBAC_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude"))$r.squared
PBACExecLong / (1-PBACExecLong)


##PBAC Total Memory Scale## 

#Cross Sectional Late LEQ
summary(lm(LEQ_t1$TotalMemoryScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$DD_PBAC_YR, na.action = "na.exclude"))
ggplot(data=LEQ_t1, aes(y=LEQ_t1$TotalMemoryScale, x=LEQ_t1$Late_LEQ, color=Group)) + geom_line(aes(group=INDDID), linetype=3, size=1) + geom_point(aes(size=1)) + guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2) + theme_classic(base_size = 17)
lm.beta::lm.beta(lm(LEQ_t1$TotalMemoryScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$DD_PBAC_YR, na.action = "na.exclude"))
plot(lm(LEQ_t1$TotalMemoryScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$DD_PBAC_YR, na.action = "na.exclude"))
shapiro.test((lm(LEQ_t1$TotalMemoryScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$DD_PBAC_YR, na.action = "na.exclude")$residuals))
PBACMemCS <- summary(lm(LEQ_t1$TotalMemoryScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$DD_PBAC_YR, na.action = "na.exclude"))$r.squared
PBACMemCS / (1-PBACMemCS)

#Longitudinal
ggplot(data=LEQ_Long, aes(y=TotalMemoryScale, x=DD_PBAC_YR)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)

PBAC_DTslope<-lmer(TotalMemoryScale~(DD_PBAC_YR|INDDID)+ SexFormatted + AgeatTestPBAC_base, data=LEQ_Long, na.action = "na.exclude", control=lmerControl(check.nobs.vs.nRE="ignore")) 
summary(PBAC_DTslope)
PBAC_DTslopes<-coef(PBAC_DTslope)$INDDID
PBAC_DTslopes$DD_MRI
PBAC_DTslopes<-setNames(cbind(rownames(PBAC_DTslopes),PBAC_DTslopes,row.names=NULL),c("INDDID", "PBAC_DTSlopes","PBAC_DTIntercept","NA1","NA2"))
PBAC_DTkeep<-PBAC_DTslopes
LEQ_Long_slopes<-merge(LEQ_Long_t1, PBAC_DTkeep, by="INDDID",all.x = TRUE)

##Late##
summary(lm(LEQ_Long_slopes$PBAC_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude"))
ggplot(data=LEQ_Long_slopes, aes(y=PBAC_DTSlopes, x=LEQ_Long_t1$Late_LEQ)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20) + labs(x="Late-Life LEQ", y="Slope of PBACC Descriptor Score")
lm.beta::lm.beta(lm(LEQ_Long_slopes$PBAC_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude"))
plot(lm(LEQ_Long_slopes$PBAC_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude"))
shapiro.test(lm(LEQ_Long_slopes$PBAC_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude")$residuals)
PBACMemLong <- summary(lm(LEQ_Long_slopes$PBAC_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude"))$r.squared
PBACMemLong / (1-PBACMemLong)
