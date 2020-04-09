library(readxl)
library(ggplot2)
library(lmerTest)
library(dplyr)
library(tidyr)
library(PMCMRplus)
library(repolr)
library(lm.beta)


## SETUP ##

setwd("/Users/nikol/Documents/LEQ PRoject/Paper 2-2020/Corey-Edits-Followup/Cog Stuff/")
LEQ<-read_xlsx("./LEQ_Data_10.1.19(3).xlsx", sheet="Sheet1")
LEQ$INDDID<-as.factor(LEQ$INDDID)
LEQ$Group<-as.factor(LEQ$Group)
LEQ$Sex<-as.factor(LEQ$Sex)
LEQ$Visit<-as.factor(LEQ$Visit)
LEQ$Long<-as.factor(LEQ$Long)
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

##Benson##

##Cross Sectional Copy+Recall Late LEQ##
summary(lm(LEQ_t1$BensonCopy~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_Benson + LEQ_t1$Overall_DD_MRI, na.action="na.exclude"))
lm.beta::lm.beta(lm(LEQ_t1$BensonCopy~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_Benson + LEQ_t1$Overall_DD_MRI, na.action="na.exclude"))

summary(lm(LEQ_t1$BensonRecall~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_Benson + LEQ_t1$Overall_DD_MRI, na.action="na.exclude"))
lm.beta::lm.beta(lm(LEQ_t1$BensonRecall~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_Benson + LEQ_t1$Overall_DD_MRI, na.action="na.exclude"))



#Longitudinal
ggplot(data=LEQ_Long, aes(y=BensonCopy, x=Overall_DD_MRI)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)+labs(x="Disease Duration (years)", y="Benson Copy")
ggplot(data=LEQ_Long, aes(y=BensonRecall, x=Overall_DD_MRI)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)+labs(x="Disease Duration (years)", y="Benson Recall")

##Copy##
BensonCopyslope<-lmer(BensonCopy~Overall_DD_MRI + (DD_MRI|INDDID)+ AgeatTest_Benson+ SexFormatted, data=LEQ_Long, na.action = "na.exclude", control=lmerControl(check.nobs.vs.nRE="ignore"), REML = FALSE) 
summary(BensonCopyslope)
BensonCopyslopes<-coef(BensonCopyslope)$INDDID
BensonCopyslopes$DD_MRI
BensonCopyslopes<-setNames(cbind(rownames(BensonCopyslopes),BensonCopyslopes,row.names=NULL),c("INDDID", "BensonCopySlopes","BensonCopyIntercept","NA1","NA2","NA3"))
BensonCopykeep<-select (BensonCopyslopes,-c(NA1,NA2,NA3))
LEQ_Long_slopes<-merge(LEQ_Long_t1, BensonCopykeep, by="INDDID",all.x = TRUE)

##Late##
summary(lm(LEQ_Long_slopes$BensonCopySlopes~LEQ_Long_t1$Late_LEQ))
lm.beta::lm.beta(lm(LEQ_Long_slopes$BensonCopySlopes~LEQ_Long_t1$Late_LEQ))



##Recall##
BensonRecallslope<-lmer(BensonRecall~Overall_DD_MRI + (DD_MRI|INDDID)+ AgeatTest_Benson+ SexFormatted, data=LEQ_Long, na.action = "na.exclude", control=lmerControl(check.nobs.vs.nRE="ignore")) 
summary(BensonRecallslope)
BensonRecallslopes<-coef(BensonRecallslope)$INDDID
BensonRecallslopes$DD_MRI
BensonRecallslopes<-setNames(cbind(rownames(BensonRecallslopes),BensonRecallslopes,row.names=NULL),c("INDDID", "BensonRecallSlopes","BensonRecallIntercept","NA1","NA2","NA3"))
BensonRecallkeep<-select (BensonRecallslopes,-c(NA1,NA2,NA3))
LEQ_Long_slopes<-merge(LEQ_Long_t1, BensonRecallkeep, by="INDDID",all.x = TRUE)

#Late#
summary(lm(LEQ_Long_slopes$BensonRecallSlopes~LEQ_Long_t1$Late_LEQ))
lm.beta::lm.beta(lm(LEQ_Long_slopes$BensonRecallSlopes~LEQ_Long_t1$Late_LEQ))



## MMSE ##

#Cross Sectional Late LEQ
summary(lm(LEQ_t1$MMSETotal~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_MMSE + LEQ_t1$Overall_DD_MRI))
ggplot(data=LEQ_t1, aes(y=LEQ_t1$MMSETotal, x=LEQ_t1$Late_LEQ, color=Group)) + geom_line(aes(group=INDDID), linetype=3, size=1) + geom_point(aes(size=1)) + guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2) + theme_classic(base_size = 17)
lm.beta::lm.beta(lm(LEQ_t1$MMSETotal~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_MMSE + LEQ_t1$Overall_DD_MRI))



#Longitudinal
ggplot(data=LEQ_Long, aes(y=MMSETotal, x=Overall_DD_MRI)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)+labs(x="Disease Duration (years)", y="MMSE")

MMSEslope<-lmer(MMSETotal~Overall_DD_MRI + (DD_MRI|INDDID)+ AgeatTest_MMSE + SexFormatted, data=LEQ_Long, control=lmerControl(check.nobs.vs.nRE="ignore")) 
summary(MMSEslope)
MMSEslopes<-coef(MMSEslope)$INDDID
MMSEslopes$DD_MRI
MMSEslopes<-setNames(cbind(rownames(MMSEslopes),MMSEslopes,row.names=NULL),c("INDDID", "MMSESlopes","MMSEIntercept","NA1","NA2","NA3"))
MMSEkeep<-select (MMSEslopes,-c(NA1,NA2,NA3))
LEQ_Long_slopes<-merge(LEQ_Long_t1, MMSEkeep, by="INDDID",all.x = TRUE)

##Late##
summary(lm(LEQ_Long_slopes$MMSESlopes~LEQ_Long_t1$Late_LEQ))
lm.beta::lm.beta(lm(LEQ_Long_slopes$MMSESlopes~LEQ_Long_t1$Late_LEQ))



## Social Behavior Observer (SBO) Checklist - Descriptor Total ##

#Cross Sectional Late LEQ
summary(lm(LEQ_t1$DescriptorTotal~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_SBO + LEQ_t1$Overall_DD_MRI, na.action = "na.exclude"))
ggplot(data=LEQ_t1, aes(y=LEQ_t1$DescriptorTotal, x=LEQ_t1$Late_LEQ, color=Group)) + geom_line(aes(group=INDDID), linetype=3, size=1) + geom_point(aes(size=1)) + guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2) + theme_classic(base_size = 17)
lm.beta::lm.beta(lm(LEQ_t1$DescriptorTotal~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_SBO + LEQ_t1$Overall_DD_MRI, na.action = "na.exclude"))


#Longitudinal
ggplot(data=LEQ_Long, aes(y=DescriptorTotal, x=Overall_DD_MRI)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)

SBO_DTslope<-lmer(DescriptorTotal~Overall_DD_MRI + (DD_MRI|INDDID)+ AgeatTest_SBO + SexFormatted, data=LEQ_Long, na.action = "na.exclude", control=lmerControl(check.nobs.vs.nRE="ignore")) 
summary(SBO_DTslope)
SBO_DTslopes<-coef(SBO_DTslope)$INDDID
SBO_DTslopes$DD_MRI
SBO_DTslopes<-setNames(cbind(rownames(SBO_DTslopes),SBO_DTslopes,row.names=NULL),c("INDDID", "SBO_DTSlopes","SBO_DTIntercept","NA1","NA2","NA3"))
SBO_DTkeep<-select (SBO_DTslopes,-c(NA1,NA2,NA3))
LEQ_Long_slopes<-merge(LEQ_Long_t1, SBO_DTkeep, by="INDDID",all.x = TRUE)

##Late##
summary(lm(LEQ_Long_slopes$SBO_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude"))
ggplot(data=LEQ_Long_slopes, aes(y=SBO_DTSlopes, x=LEQ_Long_t1$Late_LEQ)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20) + labs(x="Late-Life LEQ", y="Slope of SBOC Descriptor Score")
lm.beta::lm.beta(lm(LEQ_Long_slopes$SBO_DTSlopes~LEQ_Long_t1$Late_LEQ, na.action = "na.exclude"))



## PBAC Behavior ##

#Cross Sectional Late LEQ
summary(lm(LEQ_t1$BehavioralScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$Overall_DD_MRI))
ggplot(data=LEQ_t1, aes(y=LEQ_t1$BehavioralScale, x=LEQ_t1$Late_LEQ)) + geom_line(aes(group=INDDID), linetype=3, size=1) + geom_point(aes(size=1)) + guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2) + theme_classic(base_size = 17) + labs(x="Late-Life LEQ", y="PBAC Behavior Score")
lm.beta::lm.beta(lm(LEQ_t1$BehavioralScale~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_PBAC + LEQ_t1$Overall_DD_MRI))



#Longitudinal
ggplot(data=LEQ_Long, aes(y=BehavioralScale, x=Overall_DD_MRI)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)

PBACBehslope<-lmer(BehavioralScale~Overall_DD_MRI + (DD_MRI|INDDID)+ AgeatTest_PBAC + SexFormatted, data=LEQ_Long, control=lmerControl(check.nobs.vs.nRE="ignore")) 
summary(PBACBehslope)
PBACBehslopes<-coef(PBACBehslope)$INDDID
PBACBehslopes$DD_MRI
PBACBehslopes<-setNames(cbind(rownames(PBACBehslopes),PBACBehslopes,row.names=NULL),c("INDDID", "PBACBehslopes","PBACBehIntercept","NA1","NA2","NA3"))
PBACBehkeep<-select (PBACBehslopes,-c(NA1,NA2,NA3))
LEQ_Long_slopes<-merge(LEQ_Long_t1, PBACBehkeep, by="INDDID",all.x = TRUE)

##Late##
summary(lm(LEQ_Long_slopes$PBACBehslopes~LEQ_Long_t1$Late_LEQ))
ggplot(data=LEQ_Long_slopes, aes(y=PBACBehslopes, x=LEQ_Long_t1$Late_LEQ)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)
lm.beta::lm.beta(lm(LEQ_Long_slopes$PBACBehslopes~LEQ_Long_t1$Late_LEQ))



## BNT ##

#Cross Sectional Late LEQ
summary(lm(LEQ_t1$BNT~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_BNT + LEQ_t1$Overall_DD_MRI))
lm.beta::lm.beta(lm(LEQ_t1$BNT~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_BNT + LEQ_t1$Overall_DD_MRI))


#Longitudinal
ggplot(data=LEQ_Long, aes(y=BNT, x=Overall_DD_MRI)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)+labs(x="Disease Duration (years)", y="BNT")

BNTslope<-lmer(BNT~Overall_DD_MRI + (DD_MRI|INDDID)+ AgeatTest_BNT+ SexFormatted, data=LEQ_Long, na.action = "na.exclude", control=lmerControl(check.nobs.vs.nRE="ignore")) 
summary(BNTslope)
BNTslopes<-coef(BNTslope)$INDDID
BNTslopes$DD_MRI
BNTslopes<-setNames(cbind(rownames(BNTslopes),BNTslopes,row.names=NULL),c("INDDID", "BNTslopes","BNTIntercept","NA1","NA2","NA3"))
BNTkeep<-select (BNTslopes,-c(NA1,NA2,NA3))
LEQ_Long_slopes<-merge(LEQ_Long_t1, BNTkeep, by="INDDID",all.x = TRUE)

##Late##
summary(lm(LEQ_Long_slopes$BNTslopes~LEQ_Long_t1$Late_LEQ))
lm.beta::lm.beta(lm(LEQ_Long_slopes$BNTslopes~LEQ_Long_t1$Late_LEQ))



## F Words ##

#Cross Sectional Late LEQ
summary(lm(LEQ_t1$FWordsCorrect~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_VF + LEQ_t1$Overall_DD_MRI))
lm.beta::lm.beta(lm(LEQ_t1$FWordsCorrect~LEQ_t1$Late_LEQ + LEQ_t1$SexFormatted + LEQ_t1$AgeatTest_VF + LEQ_t1$Overall_DD_MRI))




#Longitudinal
ggplot(data=LEQ_Long, aes(y=FWordsCorrect, x=Overall_DD_MRI)) + geom_line(aes(group=INDDID), linetype=3) + geom_point(aes(size=1))+guides(size=FALSE) + geom_smooth(method = "lm", se=FALSE, size=2.5)+theme_classic(base_size = 20)+labs(x="Disease Duration (years)", y="F Words")

Fslope<-lmer(FWordsCorrect~Overall_DD_MRI + (DD_MRI|INDDID)+ AgeatTest_VF+ SexFormatted, data=LEQ_Long, na.action = "na.exclude", control=lmerControl(check.nobs.vs.nRE="ignore"))
summary(Fslope)
Fslopes<-coef(Fslope)$INDDID
Fslopes$DD_MRI
Fslopes<-setNames(cbind(rownames(Fslopes),Fslopes,row.names=NULL),c("INDDID", "Fslopes","FIntercept","NA1","NA2","NA3"))
Fkeep<-select (Fslopes,-c(NA1,NA2,NA3))
LEQ_Long_slopes<-merge(LEQ_Long_t1, Fkeep, by="INDDID",all.x = TRUE)

#Late#
summary(lm(LEQ_Long_slopes$Fslopes~LEQ_Long_t1$Late_LEQ))
lm.beta::lm.beta(lm(LEQ_Long_slopes$Fslopes~LEQ_Long_t1$Late_LEQ))








