library(ANTsR)
library(lme4)
library(reshape2)
library(lmerTest)
library(MatchIt)
library(readxl)

setwd("/Users/nikol/Documents/LEQ PRoject/Paper 2-2020/May2020/Summary_Lauren_5-15/")

##loading data, format data and generate label numbers from QuANTs output##
load("/Users/nikol/Documents/LEQ PRoject/Paper 2-2020/May2020/LS125thick.Rdata")
bvFTDLong.df <- read_excel("/Users/nikol/Documents/LEQ PRoject/Paper 2-2020/May2020/Final_LEQ_Cohort.xlsx",sheet="bvFTD-Long")

##Check QuANTs output for NAs and correct dimensions##
dim(CSLS125.df)
dim(LongLS125.df)
dim(BVvsCTRLLS125.df)
sum(is.na(CSLS125.df))
sum(is.na(LongLS125.df))
sum(is.na(BVvsCTRLLS125.df))

##check INDD demographic data predictor variables for NAs##
sum(is.na(bvFTDCS.df$AgeatMRI))
sum(is.na(bvFTDCS.df$DxDurYRS))
sum(is.na(bvFTDCS.df$SexFormatted))
sum(is.na(bvFTDCS.df$Late.Life))

sum(is.na(bvFTDLong.df$BaseDxDurYR))
sum(is.na(bvFTDLong.df$TimeYR))

##merge lausanne125 label values and demographic data from INDD##
CSLS125.df <- merge(bvFTDCS.df, CSLS125.df, by=c("id","date"))
BVvsCTRLLS125.df <- merge(BVvsCTRLLS125.df, BVvsCTRL.df, by=c("id","date"))
LongLS125.df <- merge(bvFTDLong.df, LongLS125.df, by=c("id","date"))

##Formatting##
CSLS125.df$id <- as.factor(CSLS125.df$id)
LongLS125.df$id <- as.factor(LongLS125.df$id)
BVvsCTRLLS125.df$id <- as.factor(BVvsCTRLLS125.df$id)

##generate lausanne125 label numbers
labelsgen <- colnames(CSLS125.df)
labelsgen <- labelsgen[64:282]
labelsgen <- strsplit(labelsgen,split="_")

labels <- vector(length=length(labels), mode="character")

for(i in 1:length(labelsgen)){
  labels[i] <- labelsgen[[i]][2]
}


##generate baseline cohort from longitudinal cohort##
tp1 <- subset(LongLS125.df, LongLS125.df$TimeYR==0)

##divide cohort into quantitative groups for BV vs CTRL analysis##
#BVvsCTRLLS125.df$Group <- ifelse(grepl("Normal",BVvsCTRLLS125.df$ClinicalPhenotype1),1,-1)
BVvsCTRLLS125.df$Group<-factor(BVvsCTRLLS125.df$ClinicalPhenotype1, levels = c("Normal", "bvFTD"))
##summary statistics for table 1##
CTRL.df <- subset(BVvsCTRLLS125.df, BVvsCTRLLS125.df$ClinicalPhenotype1=="Normal")

##Table1 Statistics##
sum(CSLS125.df$SexFormatted)
sum(CTRL.df$SexFormatted)
summary(CSLS125.df$AgeatMRI)
summary(CTRL.df$AgeatMRI)
summary(CSLS125.df$DxDurYRS)
summary(CSLS125.df$Education)
summary(CTRL.df$Education)

##Merge in MMSE data for MMSE Table 1 summary - Controls##
ctrlmmse.df <- merge(CTRL.df, AllControls.df, by="id")
summary(ctrlmmse.df$MMSETotal)

##Merge in MMSE data for MMSE Table 1 summary - bvFTD##
ctrlmmse.df <- data.frame(ctrlmmse.df$id, ctrlmmse.df$MMSETotal, ctrlmmse.df$Group)
cog.df <- read_excel("./Cog_Stuff/LEQ_Data_10.1.19(3).xlsx", sheet=1)
cog_t1<-subset(cog.df, Visit=="1")
bvmmse.df <- data.frame(cog_t1$INDDID, cog_t1$MMSETotal)
bvmmse.df$Group <- "bvFTD"
summary(bvmmse.df$MMSE)
colnames(ctrlmmse.df) <- c("id","MMSE", "Group")
colnames(bvmmse.df) <- c("id","MMSE", "Group")
mmsetotal.df <- rbind(ctrlmmse.df, bvmmse.df)
summary(mmsetotal.df$MMSE)

##Wilcox Testing between bv and CTRL predictor variables##

wilcox.test(AgeatMRI~Group, data=BVvsCTRLLS125.df)
wilcox.test(Education~Group, data=BVvsCTRLLS125.df)
wilcox.test(MMSE~Group, data=mmsetotal.df)

##chi sq test for sex differences##
sexmat <- matrix(c(22,9,21,10),nrow=2,ncol=2)
chisq.test(sexmat)

##regression by effect of group...this tells significant differences between groups##

##establish empty list to fill with results of lm and empty data frame to be filled with stats##
DiffList <- vector(length=length(labels), mode="list")
Diffstats <- data.frame()

##for loop...iterates thru data frame of LS125 labels for bvFTD and Controls regressing by patient group, then puts stats into a data frame##
for(i in 1:length(labels)) {
  DiffList[[i]] <- lm(as.formula(paste('lausanne125', labels[i], 'thickness_mean~Group+AgeatMRI+SexFormatted', sep="_")), data=BVvsCTRLLS125.df)
  temp <- DiffList[[i]]
  summtemp <- summary(temp)
  coef <- coef(temp)
  summcoef <- coef(summtemp)
  tempframe <- data.frame(intercept=coef[1], GroupSlope=coef[2], GroupTstat=summcoef[10], GroupP=summcoef[14], AgeSlope=coef[3], AgeTstat=summcoef[11], AgeP=summcoef[15], SexSlope=coef[4], SexTstat=summcoef[12], SexP=summcoef[16], lausanne125label=labels[i])
  Diffstats <- rbind(Diffstats, tempframe)
}


##order by uncorrected p value##
Diffstats <- Diffstats[order(Diffstats$GroupP),]

##P-value correction and re-order##
Diffstats$FWEP <- p.adjust(Diffstats$GroupP, 'bonferroni')

##take only significant p values for output##
sigDiffstats_FWE <- subset(Diffstats, Diffstats$FWEP < 0.05)

#new visualization method Chris's script##
sigDiffstats_FWE <- data.frame(sigDiffstats_FWE$lausanne125label, sigDiffstats_FWE$GroupTstat)
write.table(sigDiffstats_FWE, file="BVvsCTRLstats_FWE_tstat.txt", col.names = FALSE, row.names = FALSE)

##Create Mask of Labels Significantly Atrophied (FWEp<0.05)##
labels_FWEmask <- as.character(sigDiffstats_FWE$sigDiffstats_FWE.lausanne125label)


##Late LEQ##

##empty list to be filled with results of lm. empty dataframe for stats##
LateLEQLMlist <- vector(length=length(labels_FWEmask), mode='list')
LateLEQLMstats <- data.frame()

##for loop that iterates thru the data frame of LS125 labels_FWEmask, LEQ scores, and covariates and creates a LM for LEQ predicting LS125 label (one for each LEQ stage)##
for(i in 1:length(labels_FWEmask)) {
  LateLEQLMlist[[i]] <- lm(as.formula(paste('lausanne125', labels_FWEmask[i], 'thickness_mean~Late.Life+AgeatMRI+SexFormatted+DxDurYRS', sep="_")), data=CSLS125.df)
  tmp <- (LateLEQLMlist[[i]])
  summtmp <- summary(tmp)
  coef <- coef(tmp)
  summcoef <- coef(summtmp)
  tmpframe <- data.frame(intercept=coef[1], LateLEQSlope=coef[2], LateLEQTstat=summcoef[12], LateLEQP=summcoef[17], AgeSlope=coef[3], AgeTstat=summcoef[13], AgeP=summcoef[18], SexSlope=coef[4], SexTstat=summcoef[14], SexP=summcoef[19], DxDurSlope=coef[5], DxDurTstat=summcoef[15], DxDurP=summcoef[20], lausanne125label=labels_FWEmask[i])
  LateLEQLMstats <- rbind(LateLEQLMstats, tmpframe)
}

##order by uncorrected p value##
LateLEQLMstats <- LateLEQLMstats[order(LateLEQLMstats$LateLEQP),]

##FDR correction##
LateLEQLMstats$FDRp<-p.adjust(LateLEQLMstats$LateLEQP,'fdr')
LateLEQLMstats$FWep <- p.adjust(LateLEQLMstats$LateLEQP, 'bonferroni')

#output significant image##
sigLatestats_FDR <- subset(LateLEQLMstats, LateLEQLMstats$FDRp<0.05)
sigLatestats_FDR <- subset(sigLatestats, sigLatestats$LateLEQSlope>0)

sigLatestats_unc005 <- subset(LateLEQLMstats, LateLEQLMstats$LateLEQP<0.005)
sigLatestats_unc005 <- subset(sigLatestats_unc005, sigLatestats_unc005$LateLEQSlope>0)

##new rendering format for Chris' Lausanne Render Script##
sigLatestats_FDR <- data.frame(sigLatestats$lausanne125label, sigLatestats$FDRp)
write.table(sigLatestats_FDR, file="LateCS_FDR.txt", col.names = FALSE, row.names = FALSE)

sigLatestats_unc005 <- data.frame(sigLatestats_unc005$lausanne125label, sigLatestats_unc005$LateLEQTstat)
write.table(sigLatestats_unc005, file="LateCS_unc005_tstat.txt", col.names = FALSE, row.names = FALSE)

##Early LEQ##

##empty list to be filled with results of lm. empty dataframe for stats##
EarlyLEQLMlist <- vector(length=length(labels_FWEmask), mode='list')
EarlyLEQLMstats <- data.frame()

##for loop that iterates thru the data frame of LS125 labels_FWEmask, LEQ scores, and covariates and creates a LM for LEQ predicting LS125 label (one for each LEQ stage)##
for(i in 1:length(labels_FWEmask)) {
  EarlyLEQLMlist[[i]] <- lm(as.formula(paste('lausanne125', labels_FWEmask[i], 'thickness_mean~Young.Adulthood+AgeatMRI+SexFormatted+DxDurYRS', sep="_")), data=CSLS125.df)
  tmp <- (EarlyLEQLMlist[[i]])
  summtmp <- summary(tmp)
  coef <- coef(tmp)
  summcoef <- coef(summtmp)
  tmpframe <- data.frame(intercept=coef[1], EarlyLEQSlope=coef[2], EarlyLEQTstat=summcoef[12], EarlyLEQP=summcoef[17], AgeSlope=coef[3], AgeTstat=summcoef[13], AgeP=summcoef[18], SexSlope=coef[4], SexTstat=summcoef[14], SexP=summcoef[19], DxDurSlope=coef[5], DxDurTstat=summcoef[15], DxDurP=summcoef[20], lausanne125label=labels_FWEmask[i])
  EarlyLEQLMstats <- rbind(EarlyLEQLMstats, tmpframe)
}

##order by uncorrected p value##
EarlyLEQLMstats <- EarlyLEQLMstats[order(EarlyLEQLMstats$EarlyLEQP),]

##FDR correction##
EarlyLEQLMstats$FDRp<-p.adjust(EarlyLEQLMstats$EarlyLEQP,'fdr')
EarlyLEQLMstats$FWep <- p.adjust(EarlyLEQLMstats$EarlyLEQP, 'bonferroni')

#output significant image##
sigEarlystats <- subset(EarlyLEQLMstats, EarlyLEQLMstats$FDRp<0.05)
sigEarlystats <- subset(sigEarlystats, sigEarlystats$EarlyLEQSlope>0)

##new rendering format for Chris' Lausanne Render Script##
sigEarlystats <- data.frame(sigEarlystats$lausanne125label, sigEarlystats$FDRp)
write.table(sigEarlystats, file="EarlyCSStats.txt", col.names = FALSE, row.names = FALSE)

##Mid LEQ##

##empty list to be filled with results of lm. empty dataframe for stats##
MidLEQLMlist <- vector(length=length(labels_FWEmask), mode='list')
MidLEQLMstats <- data.frame()

##for loop that iterates thru the data frame of LS125 labels_FWEmask, LEQ scores, and covariates and creates a LM for LEQ predicting LS125 label (one for each LEQ stage)##
for(i in 1:length(labels_FWEmask)) {
  MidLEQLMlist[[i]] <- lm(as.formula(paste('lausanne125', labels_FWEmask[i], 'thickness_mean~Mid.Life+AgeatMRI+SexFormatted+DxDurYRS', sep="_")), data=CSLS125.df)
  tmp <- (MidLEQLMlist[[i]])
  summtmp <- summary(tmp)
  coef <- coef(tmp)
  summcoef <- coef(summtmp)
  tmpframe <- data.frame(intercept=coef[1], MidLEQSlope=coef[2], MidLEQTstat=summcoef[12], MidLEQP=summcoef[17], AgeSlope=coef[3], AgeTstat=summcoef[13], AgeP=summcoef[18], SexSlope=coef[4], SexTstat=summcoef[14], SexP=summcoef[19], DxDurSlope=coef[5], DxDurTstat=summcoef[15], DxDurP=summcoef[20], lausanne125label=labels_FWEmask[i])
  MidLEQLMstats <- rbind(MidLEQLMstats, tmpframe)
}

##order by uncorrected p value##
MidLEQLMstats <- MidLEQLMstats[order(MidLEQLMstats$MidLEQP),]

##FDR correction##
MidLEQLMstats$FDRp<-p.adjust(MidLEQLMstats$MidLEQP,'fdr')
MidLEQLMstats$FWep <- p.adjust(MidLEQLMstats$MidLEQP, 'bonferroni')

#output significant image##
sigMidstats <- subset(MidLEQLMstats, MidLEQLMstats$FDRp<0.05)
sigMidstats <- subset(sigMidstats, sigMidstats$MidLEQSlope>0)

##new rendering format for Chris' Lausanne Render Script##
sigMidstats <- data.frame(sigMidstats$lausanne125label, sigMidstats$FDRp)
write.table(sigMidstats, file="MidCSStats.txt", col.names = FALSE, row.names = FALSE)


##Longitudinal##

##empty list to be filled with results of lmer. empty dataframe for slopes##
LMERlist <- vector(length=length(labels_FWEmask), mode='list')
LMERslopes <- data.frame()
singtest <- data.frame()

##for loop generating model for each LS125 label##
for(i in 1:length(labels_FWEmask)) {
  LMERlist[[i]] <- lmer(as.formula(paste('lausanne125', labels_FWEmask[i], 'thickness_mean~(DxDurYRS|id)+AgeatMRI_base+SexFormatted', sep="_")), data=LongLS125.df, na.action = "na.exclude")
  sing <- isSingular(LMERlist[[i]])
  tempsing <- data.frame(label=paste('lausanne125',labels_FWEmask[i],sep="_"), sing=sing)
  singtest <- rbind(singtest,tempsing)
  tmp <- (LMERlist[[i]])
  placeholder <- coef(tmp)
  temp <- placeholder$id
  slopes <- setNames(cbind(rownames(temp), temp,row.names=NULL),c("id","Slopes","Intercept","NA1","NA2"))
  tmpframe <- data.frame(id=slopes$id, slope=slopes$Slopes, lausanne125label=paste('lausanne125',labels_FWEmask[i],'slope',sep="_"))
  LMERslopes <- rbind(LMERslopes, tmpframe)
}

##reformat for further testing and add demographic values from INDD##
SlopesforTesting <- dcast(LMERslopes, id~lausanne125label, value.var="slope")
SlopesforTesting <- merge(SlopesforTesting, tp1, by="id")

##Relate SLopes to LEQ -- Late##
LateLongCorList <- vector(length=length(labels_FWEmask), mode='list')
LateLongLMstats <- data.frame()

##relate all LS125 region rate of change to late life LEQ##
for(i in 1:length(labels_FWEmask)){
  LateLongCorList[[i]] <- lm(as.formula(paste('lausanne125', labels_FWEmask[i], 'slope~Late.Life', sep="_")), data=SlopesforTesting)
  tmp <- (LateLongCorList[[i]])
  summcoef <- coef(summary(tmp))
  tmpframe <- data.frame(LateLEQP=summcoef[8], LateLEQSlope=summcoef[2], LateLEQTstat=summcoef[6],lausanne125label=labels_FWEmask[i])
  LateLongLMstats <- rbind(LateLongLMstats, tmpframe)
}

##order dataframe and p-value correction##
LateLongLMstats <- LateLongLMstats[order(LateLongLMstats$LateLEQP),]
LateLongLMstats$FDRp <- p.adjust(LateLongLMstats$LateLEQP, 'fdr')
LateLongLMstats$FWEp <- p.adjust(LateLongLMstats$LateLEQP, 'bonferroni')

unc005_LateLongLMstats <- subset(LateLongLMstats, LateLongLMstats$LateLEQP<0.005)
unc005_LateLongLMstats <- data.frame(unc005_LateLongLMstats$lausanne125label, unc005_LateLongLMstats$LateLEQTstat)
write.table(unc005_LateLongLMstats, file="unc005_Coreymodel_tstat.txt", col.names = FALSE, row.names = FALSE)  

##Mid-Life##
MidLongCorList <- vector(length=length(labels_FWEmask), mode='list')
MidLongLMstats <- data.frame()

##reMid all LS125 region rate of change to Mid life LEQ##
for(i in 1:length(labels_FWEmask)){
  MidLongCorList[[i]] <- lm(as.formula(paste('lausanne125', labels_FWEmask[i], 'slope~Mid.Life', sep="_")), data=SlopesforTesting)
  tmp <- (MidLongCorList[[i]])
  summcoef <- coef(summary(tmp))
  tmpframe <- data.frame(MidLEQP=summcoef[8], MidLEQSlope=summcoef[2], lausanne125label=labels_FWEmask[i])
  MidLongLMstats <- rbind(MidLongLMstats, tmpframe)
}

##order dataframe and p-value correction##
MidLongLMstats <- MidLongLMstats[order(MidLongLMstats$MidLEQP),]
MidLongLMstats$FDRp <- p.adjust(MidLongLMstats$MidLEQP, 'fdr')
MidLongLMstats$FWEp <- p.adjust(MidLongLMstats$MidLEQP, 'bonferroni')

##Early-Life##
EarlyLongCorList <- vector(length=length(labels_FWEmask), mode='list')
EarlyLongLMstats <- data.frame()

##reEarly all LS125 region rate of change to Early life LEQ##
for(i in 1:length(labels_FWEmask)){
  EarlyLongCorList[[i]] <- lm(as.formula(paste('lausanne125', labels_FWEmask[i], 'slope~Young.Adulthood', sep="_")), data=SlopesforTesting)
  tmp <- (EarlyLongCorList[[i]])
  summcoef <- coef(summary(tmp))
  tmpframe <- data.frame(EarlyLEQP=summcoef[8], EarlyLEQSlope=summcoef[2], lausanne125label=labels_FWEmask[i])
  EarlyLongLMstats <- rbind(EarlyLongLMstats, tmpframe)
}

##order dataframe and p-value correction##
EarlyLongLMstats <- EarlyLongLMstats[order(EarlyLongLMstats$EarlyLEQP),]
EarlyLongLMstats$FDRp <- p.adjust(EarlyLongLMstats$EarlyLEQP, 'fdr')
EarlyLongLMstats$FWEp <- p.adjust(EarlyLongLMstats$EarlyLEQP, 'bonferroni')
