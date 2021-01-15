library(lme4)
library(reshape2)
library(lmerTest)
library(MatchIt)
library(readxl)
library(corrplot)
library(Hmisc)
library(lm.beta)
library(yhat)

setwd("/Users/nikol/Documents/LEQ PRoject/Paper 2-2020/May2020/ReviewerComments_11-18/")

##loading data, format data and generate label numbers from QuANTs output##
load("/Users/nikol/Documents/LEQ PRoject/Paper 2-2020/May2020/LS125thick.Rdata")
bvFTDLong.df <- read_excel("/Users/nikol/Documents/LEQ PRoject/Paper 2-2020/May2020/ReviewerComments_11-18/Final_LEQ_Cohort.xlsx",sheet="bvFTD-Long")

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
cog.df <- read_excel("/Users/nikol/Documents/LEQ PRoject/Paper 2-2020/May2020/ReviewerComments_11-18/Cog_Stuff/LEQ_Data_10.1.19(3).xlsx", sheet=1)
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

##Correlation Matrix for LEQ vars + Education##
cormat.df <- data.frame(bvFTDCS.df$Young.Adulthood, bvFTDCS.df$Mid.Life, bvFTDCS.df$Late.Life, bvFTDCS.df$Education)
colnames(cormat.df) <- c("Early Life", "Mid Life", "Late Life", "Education")
cor <- cor(cormat.df)
corp <- rcorr(as.matrix(cormat.df))
corrplot(cor)
corrplot(corp$P)


##read in sequences for BVvsCTRL##
BVvsCTRLseq.df <- read.csv("BVvsCTRLseqsearch.csv")
BVvsCTRLLS125.df <- merge(BVvsCTRLLS125.df, BVvsCTRLseq.df, by=c("id","date"))

##regression by effect of group...this tells significant differences between groups##

##establish empty list to fill with results of lm and empty data frame to be filled with stats##
DiffList <- vector(length=length(labels), mode="list")
Diffstats <- data.frame()

##for loop...iterates thru data frame of LS125 labels for bvFTD and Controls regressing by patient group, then puts stats into a data frame##
##CHANGE -- removed age+sex here##
for(i in 1:length(labels)) {
  DiffList[[i]] <- lm(as.formula(paste('lausanne125', labels[i], 'thickness_mean~Group+SequenceFormatted', sep="_")), data=BVvsCTRLLS125.df)
  temp <- DiffList[[i]]
  summtemp <- summary(temp)
  coef <- coef(temp)
  summcoef <- coef(summtemp)
  beta <- lm.beta(DiffList[[i]])$standardized.coefficients[2][[1]]
  tempframe <- data.frame(lausanne125label=labels[i], Intercept=summcoef[1], GroupSlope=summcoef[2], GroupTstat=summcoef[8], GroupP=summcoef[11], GroupBeta=beta)
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
write.table(sigDiffstats_FWE, file="BVvsCTRLstats_FWE_tstat_SequenceCovary.txt", col.names = FALSE, row.names = FALSE)

##Create Mask of Labels Significantly Atrophied (FWEp<0.05)##
labels_FWEmask <- as.character(sigDiffstats_FWE$sigDiffstats_FWE.lausanne125label)




# Full Analysis w/Scanner Covariate ---------------------------------------

##Late LEQ##
##empty list to be filled with results of lm. empty dataframe for stats##
LateLEQLMlist <- vector(length=length(labels_FWEmask), mode='list')
LateLEQLMstats <- data.frame()

##for loop that iterates thru the data frame of LS125 labels_FWEmask, LEQ scores, and covariates and creates a LM for LEQ predicting LS125 label (one for each LEQ stage)##
for(i in 1:length(labels_FWEmask)) {
  LateLEQLMlist[[i]] <- lm(as.formula(paste('lausanne125', labels_FWEmask[i], 'thickness_mean~Late.Life+AgeatMRI+SexFormatted+DxDurYRS+T1ProtocolFormatted', sep="_")), data=CSLS125.df)
  tmp <- (LateLEQLMlist[[i]])
  summtmp <- summary(tmp)
  coef <- coef(tmp)
  summcoef <- coef(summtmp)
  beta <- lm.beta(LateLEQLMlist[[i]])$standardized.coefficients[2][[1]]
  adjrsq <- summtmp$adj.r.squared
  rq <- summtmp$r.squared
  cohenfsq <- rq / (1-rq)
  tmpframe <- data.frame(lausanne125label=labels_FWEmask[i], Intercept=summcoef[1], LateLEQSlope=summcoef[2], LateLEQTstat=summcoef[14], LateLEQP=summcoef[20],LateLEQBeta=beta, Rsquared=rq, fsquared = cohenfsq)
  LateLEQLMstats <- rbind(LateLEQLMstats, tmpframe)
}

##order by uncorrected p value##
LateLEQLMstats <- LateLEQLMstats[order(LateLEQLMstats$LateLEQP),]

##FDR correction##
LateLEQLMstats$FDRp<-p.adjust(LateLEQLMstats$LateLEQP,'fdr')
LateLEQLMstats$FWep <- p.adjust(LateLEQLMstats$LateLEQP, 'bonferroni')

##read in ls125 names##
ls125_names.df <- read.csv("/Users/nikol/Documents/LEQ PRoject/Paper 2-2020/May2020/ReviewerComments_11-18/Imaging_Results/Lausanne_Scale125.csv")
colnames(ls125_names.df) <- c("lausanne125label","Label.Name")

##merge ls125 names into table##
LateLEQLMstats <- merge(ls125_names.df, LateLEQLMstats, by="lausanne125label")
LateLEQLMstats <- LateLEQLMstats[order(LateLEQLMstats$LateLEQP),]

#output significant image##
sigLatestats_FWE <- subset(LateLEQLMstats, LateLEQLMstats$FWep<0.05)
sigLatestats_FWE <- subset(sigLatestats_FWE, sigLatestats_FWE$LateLEQSlope>0)

##new rendering format for Chris' Lausanne Render Script##
sigLatestats_FWE <- data.frame(sigLatestats_FWE$lausanne125label, sigLatestats_FWE$LateLEQTstat)
write.table(sigLatestats_FWE, file="LateCS_FWE_SeqCovary.txt", col.names = FALSE, row.names = FALSE)

##write out f^2 heatmap table##
LateLEQCS_f2 <- data.frame(LateLEQLMstats$lausanne125label, LateLEQLMstats$fsquared)
write.table(LateLEQCS_f2, file="LateLEQCS_F2.txt", row.names = FALSE, col.names = FALSE)

##Early LEQ##

##empty list to be filled with results of lm. empty dataframe for stats##
EarlyLEQLMlist <- vector(length=length(labels_FWEmask), mode='list')
EarlyLEQLMstats <- data.frame()

##for loop that iterates thru the data frame of LS125 labels_FWEmask, LEQ scores, and covariates and creates a LM for LEQ predicting LS125 label (one for each LEQ stage)##
for(i in 1:length(labels_FWEmask)) {
  EarlyLEQLMlist[[i]] <- lm(as.formula(paste('lausanne125', labels_FWEmask[i], 'thickness_mean~Young.Adulthood+AgeatMRI+SexFormatted+DxDurYRS+T1ProtocolFormatted', sep="_")), data=CSLS125.df)
  tmp <- (EarlyLEQLMlist[[i]])
  summtmp <- summary(tmp)
  coef <- coef(tmp)
  summcoef <- coef(summtmp)
  beta <- lm.beta(EarlyLEQLMlist[[i]])$standardized.coefficients[2][[1]]
  adjrsq <- summtmp$adj.r.squared
  tmpframe <- data.frame(lausanne125label=labels_FWEmask[i], Intercept=summcoef[1], EarlyLEQSlope=summcoef[2], EarlyLEQTstat=summcoef[14], EarlyLEQP=summcoef[20],EarlyLEQBeta=beta, Adj_Rsquared=adjrsq)
  EarlyLEQLMstats <- rbind(EarlyLEQLMstats, tmpframe)
}

##order by uncorrected p value##
EarlyLEQLMstats <- EarlyLEQLMstats[order(EarlyLEQLMstats$EarlyLEQP),]

##FDR correction##
EarlyLEQLMstats$FDRp<-p.adjust(EarlyLEQLMstats$EarlyLEQP,'fdr')
EarlyLEQLMstats$FWep <- p.adjust(EarlyLEQLMstats$EarlyLEQP, 'bonferroni')

##merge ls125 names into table##
EarlyLEQLMstats <- merge(ls125_names.df, EarlyLEQLMstats, by="lausanne125label")
EarlyLEQLMstats <- EarlyLEQLMstats[order(EarlyLEQLMstats$EarlyLEQP),]

##Mid LEQ##

##empty list to be filled with results of lm. empty dataframe for stats##
MidLEQLMlist <- vector(length=length(labels_FWEmask), mode='list')
MidLEQLMstats <- data.frame()

##for loop that iterates thru the data frame of LS125 labels_FWEmask, LEQ scores, and covariates and creates a LM for LEQ predicting LS125 label (one for each LEQ stage)##
for(i in 1:length(labels_FWEmask)) {
  MidLEQLMlist[[i]] <- lm(as.formula(paste('lausanne125', labels_FWEmask[i], 'thickness_mean~Mid.Life+AgeatMRI+SexFormatted+DxDurYRS+T1ProtocolFormatted', sep="_")), data=CSLS125.df)
  tmp <- (MidLEQLMlist[[i]])
  summtmp <- summary(tmp)
  coef <- coef(tmp)
  summcoef <- coef(summtmp)
  beta <- lm.beta(MidLEQLMlist[[i]])$standardized.coefficients[2][[1]]
  adjrsq <- summtmp$adj.r.squared
  tmpframe <- data.frame(lausanne125label=labels_FWEmask[i], Intercept=summcoef[1], MidLEQSlope=summcoef[2], MidLEQTstat=summcoef[14], MidLEQP=summcoef[20],MidLEQBeta=beta, Adj_Rsquared=adjrsq)
  MidLEQLMstats <- rbind(MidLEQLMstats, tmpframe)
}

##order by uncorrected p value##
MidLEQLMstats <- MidLEQLMstats[order(MidLEQLMstats$MidLEQP),]

##FDR correction##
MidLEQLMstats$FDRp<-p.adjust(MidLEQLMstats$MidLEQP,'fdr')
MidLEQLMstats$FWep <- p.adjust(MidLEQLMstats$MidLEQP, 'bonferroni')

##merge ls125 names into table##
MidLEQLMstats <- merge(ls125_names.df, MidLEQLMstats, by="lausanne125label")
MidLEQLMstats <- MidLEQLMstats[order(MidLEQLMstats$MidLEQP),]

##Total LEQ##

##empty list to be filled with results of lm. empty dataframe for stats##
TotalLEQLMlist <- vector(length=length(labels_FWEmask), mode='list')
TotalLEQLMstats <- data.frame()

##for loop that iterates thru the data frame of LS125 labels_FWEmask, LEQ scores, and covariates and creates a LM for LEQ predicting LS125 label (one for each LEQ stage)##
for(i in 1:length(labels_FWEmask)) {
  TotalLEQLMlist[[i]] <- lm(as.formula(paste('lausanne125', labels_FWEmask[i], 'thickness_mean~Total.LEQ+AgeatMRI+SexFormatted+DxDurYRS+T1ProtocolFormatted', sep="_")), data=CSLS125.df)
  tmp <- (TotalLEQLMlist[[i]])
  summtmp <- summary(tmp)
  coef <- coef(tmp)
  summcoef <- coef(summtmp)
  beta <- lm.beta(TotalLEQLMlist[[i]])$standardized.coefficients[2][[1]]
  adjrsq <- summtmp$adj.r.squared
  tmpframe <- data.frame(lausanne125label=labels_FWEmask[i], Intercept=summcoef[1], TotalLEQSlope=summcoef[2], TotalLEQTstat=summcoef[14], TotalLEQP=summcoef[20],TotalLEQBeta=beta, Adj_Rsquared=adjrsq)
  TotalLEQLMstats <- rbind(TotalLEQLMstats, tmpframe)
}

##order by uncorrected p value##
TotalLEQLMstats <- TotalLEQLMstats[order(TotalLEQLMstats$TotalLEQP),]

##FDR correction##
TotalLEQLMstats$FDRp<-p.adjust(TotalLEQLMstats$TotalLEQP,'fdr')
TotalLEQLMstats$FWep <- p.adjust(TotalLEQLMstats$TotalLEQP, 'bonferroni')

##merge ls125 names into table##
TotalLEQLMstats <- merge(ls125_names.df, TotalLEQLMstats, by="lausanne125label")
TotalLEQLMstats <- TotalLEQLMstats[order(TotalLEQLMstats$TotalLEQP),]

##write out tables for supplementary material##
write.csv(LateLEQLMstats, file="./Imaging_Results/Tables/CS_LateLEQ.csv", row.names = FALSE)
write.csv(EarlyLEQLMstats, file="./Imaging_Results/Tables/CS_EarlyLEQ.csv", row.names = FALSE)
write.csv(MidLEQLMstats, file="./Imaging_Results/Tables/CS_MidLEQ.csv", row.names = FALSE)
write.csv(TotalLEQLMstats, file="./Imaging_Results/Tables/CS_TotalLEQ.csv", row.names = FALSE)



##Longitudinal##

##empty list to be filled with results of lmer. empty dataframe for slopes##
LMERlist <- vector(length=length(labels_FWEmask), mode='list')
LMERslopes <- data.frame()
singtest <- data.frame()

##for loop generating model for each LS125 label##
for(i in 1:length(labels_FWEmask)) {
  LMERlist[[i]] <- lmer(as.formula(paste('lausanne125', labels_FWEmask[i], 'thickness_mean~(DxDurYRS|id)+AgeatMRI_base+SexFormatted+T1ProtocolFormatted', sep="_")), data=LongLS125.df, na.action = "na.exclude")
  sing <- isSingular(LMERlist[[i]])
  tempsing <- data.frame(label=paste('lausanne125',labels_FWEmask[i],sep="_"), sing=sing)
  singtest <- rbind(singtest,tempsing)
  tmp <- (LMERlist[[i]])
  placeholder <- coef(tmp)
  temp <- placeholder$id
  slopes <- setNames(cbind(rownames(temp), temp,row.names=NULL),c("id","Slopes","Intercept","NA1","NA2","NA3"))
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
  beta <- lm.beta(LateLongCorList[[i]])$standardized.coefficients[2][[1]]
  adjrsq <- summary(tmp)$adj.r.squared
  rq <- summary(tmp)$r.squared
  cohenf2 <- rq / (1-rq)
  tmpframe <- data.frame(lausanne125label=labels_FWEmask[i],LateLEQP=summcoef[8], LateLEQSlope=summcoef[2], LateLEQTstat=summcoef[6],LateLEQBeta=beta,Rsquared=rq,Fsquared = cohenf2)
  LateLongLMstats <- rbind(LateLongLMstats, tmpframe)
}

##order dataframe and p-value correction##
LateLongLMstats <- LateLongLMstats[order(LateLongLMstats$LateLEQP),]
LateLongLMstats$FDRp <- p.adjust(LateLongLMstats$LateLEQP, 'fdr')
LateLongLMstats$FWEp <- p.adjust(LateLongLMstats$LateLEQP, 'bonferroni')

##merge ls125 names into table##
LateLongLMstats <- merge(ls125_names.df, LateLongLMstats, by="lausanne125label")
LateLongLMstats <- LateLongLMstats[order(LateLongLMstats$LateLEQP),]

##write out output##
FWE_LateLongLMstats <- subset(LateLongLMstats, LateLongLMstats$FWEp<0.05)
FWE_LateLongLMstats <- data.frame(FWE_LateLongLMstats$lausanne125label, FWE_LateLongLMstats$LateLEQTstat)
write.table(FWE_LateLongLMstats, file="FWE_Coreymodel_tstat_SeqCovary.txt", col.names = FALSE, row.names = FALSE) 

##write out Fsquared heatmap##
LateLEQLong_f2 <- data.frame(LateLongLMstats$lausanne125label, LateLongLMstats$Fsquared)
write.table(LateLEQLong_f2, file="LateLEQLong_f2.txt", row.names = FALSE, col.names = FALSE)

##Mid-Life##
MidLongCorList <- vector(length=length(labels_FWEmask), mode='list')
MidLongLMstats <- data.frame()

##reMid all LS125 region rate of change to Mid life LEQ##
for(i in 1:length(labels_FWEmask)){
  MidLongCorList[[i]] <- lm(as.formula(paste('lausanne125', labels_FWEmask[i], 'slope~Mid.Life', sep="_")), data=SlopesforTesting)
  tmp <- (MidLongCorList[[i]])
  summcoef <- coef(summary(tmp))
  beta <- lm.beta(MidLongCorList[[i]])$standardized.coefficients[2][[1]]
  adjrsq <- summary(tmp)$adj.r.squared
  tmpframe <- data.frame(lausanne125label=labels_FWEmask[i],MidLEQP=summcoef[8], MidLEQSlope=summcoef[2], MidLEQTstat=summcoef[6], MidLEQBeta=beta,Adj_Rsquared=adjrsq)
  MidLongLMstats <- rbind(MidLongLMstats, tmpframe)
}

##order dataframe and p-value correction##
MidLongLMstats <- MidLongLMstats[order(MidLongLMstats$MidLEQP),]
MidLongLMstats$FDRp <- p.adjust(MidLongLMstats$MidLEQP, 'fdr')
MidLongLMstats$FWEp <- p.adjust(MidLongLMstats$MidLEQP, 'bonferroni')

##merge ls125 names into table##
MidLongLMstats <- merge(ls125_names.df, MidLongLMstats, by="lausanne125label")
MidLongLMstats <- MidLongLMstats[order(MidLongLMstats$MidLEQP),]

##Early-Life##
EarlyLongCorList <- vector(length=length(labels_FWEmask), mode='list')
EarlyLongLMstats <- data.frame()

##reEarly all LS125 region rate of change to Early life LEQ##
for(i in 1:length(labels_FWEmask)){
  EarlyLongCorList[[i]] <- lm(as.formula(paste('lausanne125', labels_FWEmask[i], 'slope~Young.Adulthood', sep="_")), data=SlopesforTesting)
  tmp <- (EarlyLongCorList[[i]])
  summcoef <- coef(summary(tmp))
  beta <- lm.beta(EarlyLongCorList[[i]])$standardized.coefficients[2][[1]]
  adjrsq <- summary(tmp)$adj.r.squared
  tmpframe <- data.frame(lausanne125label=labels_FWEmask[i],EarlyLEQP=summcoef[8], EarlyLEQSlope=summcoef[2],EarlyLEQTstat=summcoef[6],EarlyLEQBeta=beta,Adj_Rsquared=adjrsq)
  EarlyLongLMstats <- rbind(EarlyLongLMstats, tmpframe)
}

##order dataframe and p-value correction##
EarlyLongLMstats <- EarlyLongLMstats[order(EarlyLongLMstats$EarlyLEQP),]
EarlyLongLMstats$FDRp <- p.adjust(EarlyLongLMstats$EarlyLEQP, 'fdr')
EarlyLongLMstats$FWEp <- p.adjust(EarlyLongLMstats$EarlyLEQP, 'bonferroni')

##merge ls125 names into table##
EarlyLongLMstats <- merge(ls125_names.df, EarlyLongLMstats, by="lausanne125label")
EarlyLongLMstats <- EarlyLongLMstats[order(EarlyLongLMstats$EarlyLEQP),]

##Total LEQ##
TotalLongCorList <- vector(length=length(labels_FWEmask), mode='list')
TotalLongLMstats <- data.frame()

##reTotal all LS125 region rate of change to Total life LEQ##
for(i in 1:length(labels_FWEmask)){
  TotalLongCorList[[i]] <- lm(as.formula(paste('lausanne125', labels_FWEmask[i], 'slope~Total.LEQ', sep="_")), data=SlopesforTesting)
  tmp <- (TotalLongCorList[[i]])
  summcoef <- coef(summary(tmp))
  beta <- lm.beta(TotalLongCorList[[i]])$standardized.coefficients[2][[1]]
  adjrsq <- summary(tmp)$adj.r.squared
  tmpframe <- data.frame(lausanne125label=labels_FWEmask[i],TotalLEQP=summcoef[8], TotalLEQSlope=summcoef[2],TotalLEQTstat=summcoef[6], TotalLEQBeta=beta,Adj_Rsquared=adjrsq)
  TotalLongLMstats <- rbind(TotalLongLMstats, tmpframe)
}

##order dataframe and p-value correction##
TotalLongLMstats <- TotalLongLMstats[order(TotalLongLMstats$TotalLEQP),]
TotalLongLMstats$FDRp <- p.adjust(TotalLongLMstats$TotalLEQP, 'fdr')
TotalLongLMstats$FWEp <- p.adjust(TotalLongLMstats$TotalLEQP, 'bonferroni')

##merge ls125 names into table##
TotalLongLMstats <- merge(ls125_names.df, TotalLongLMstats, by="lausanne125label")
TotalLongLMstats <- TotalLongLMstats[order(TotalLongLMstats$TotalLEQP),]

##write out tables for supplement##
write.csv(LateLongLMstats, file="./Imaging_Results/Tables/Long_LateLEQ.csv", row.names = FALSE)
write.csv(EarlyLongLMstats, file="./Imaging_Results/Tables/Long_EarlyLEQ.csv", row.names = FALSE)
write.csv(MidLongLMstats, file="./Imaging_Results/Tables/Long_MidLEQ.csv", row.names = FALSE)
write.csv(TotalLongLMstats, file="./Imaging_Results/Tables/Long_TotalLEQ.csv", row.names = FALSE)


# LEQ * Time Longitudinal Analysis w/Scanner Covariate----------------------------------------

##Late Life##

##Initialize dataframes##
LMERlist_inter_late <- vector(length=length(labels), mode='list')
LMERslopes_inter_late <- data.frame()

##for loop generating model for each LS125 label##
for(i in 1:length(labels_FWEmask)) {
  LMERlist_inter_late[[i]] <- lmer(as.formula(paste('lausanne125', labels[i], 'thickness_mean~(1|id)+TimeYR*Late.Life+AgeatMRI_base+DxDurYRS+SexFormatted + T1ProtocolFormatted', sep="_")), data=LongLS125.df, na.action = "na.exclude")
  tmp <- coef(summary(LMERlist_inter_late[[i]]))
  tmpframe <- data.frame(inter_p = tmp[40], inter_Tstat = tmp[32], lausanne125label=labels_FWEmask[i])
  LMERslopes_inter_late <- rbind(LMERslopes_inter_late, tmpframe)
}

##P-value correcting##
LMERslopes_inter_late <- LMERslopes_inter_late[order(LMERslopes_inter_late$inter_p),]
LMERslopes_inter_late$FDRp <- p.adjust(LMERslopes_inter_late$inter_p, 'fdr')
LMERslopes_inter_late$FWEp <- p.adjust(LMERslopes_inter_late$inter_p, 'bonferroni')

sig_LMER_inter_late_unc005 <- subset(LMERslopes_inter_late, LMERslopes_inter_late$inter_p<0.005)
sig_LMER_inter_late_unc005 <- data.frame(sig_LMER_inter_late_unc005$lausanne125label, sig_LMER_inter_late_unc005$inter_Tstat)
write.table(sig_LMER_inter_late_unc005, file="Late_Interaction_SeqCovary_unc005.txt", col.names = FALSE, row.names = FALSE)

##Mid Life##

##Initialize dataframes##
LMERlist_inter_Mid <- vector(length=length(labels), mode='list')
LMERslopes_inter_Mid <- data.frame()

##for loop generating model for each LS125 label##
for(i in 1:length(labels_FWEmask)) {
  LMERlist_inter_Mid[[i]] <- lmer(as.formula(paste('lausanne125', labels[i], 'thickness_mean~(1|id)+TimeYR*Mid.Life+AgeatMRI_base+DxDurYRS+SexFormatted+T1ProtocolFormatted', sep="_")), data=LongLS125.df, na.action = "na.exclude")
  tmp <- coef(summary(LMERlist_inter_Mid[[i]]))
  tmpframe <- data.frame(inter_p = tmp[40], inter_Tstat = tmp[32], lausanne125label=labels_FWEmask[i])
  LMERslopes_inter_Mid <- rbind(LMERslopes_inter_Mid, tmpframe)
}

##P-value correcting##
LMERslopes_inter_Mid <- LMERslopes_inter_Mid[order(LMERslopes_inter_Mid$inter_p),]
LMERslopes_inter_Mid$FDRp <- p.adjust(LMERslopes_inter_Mid$inter_p, 'fdr')
LMERslopes_inter_Mid$FWEp <- p.adjust(LMERslopes_inter_Mid$inter_p, 'bonferroni')

sig_LMER_inter_mid_unc005 <- subset(LMERslopes_inter_Mid, LMERslopes_inter_Mid$inter_p<0.005)
sig_LMER_inter_mid_unc005 <- data.frame(sig_LMER_inter_mid_unc005$lausanne125label, sig_LMER_inter_mid_unc005$inter_Tstat)
write.table(sig_LMER_inter_mid_unc005, file="Mid_Interaction_SeqCovary_unc005.txt", col.names = FALSE, row.names = FALSE)

##Early Life##

##Initialize dataframes##
LMERlist_inter_Early <- vector(length=length(labels), mode='list')
LMERslopes_inter_Early <- data.frame()

##for loop generating model for each LS125 label##
for(i in 1:length(labels_FWEmask)) {
  LMERlist_inter_Early[[i]] <- lmer(as.formula(paste('lausanne125', labels[i], 'thickness_mean~(1|id)+TimeYR*Young.Adulthood+AgeatMRI_base+DxDurYRS+SexFormatted+T1ProtocolFormatted', sep="_")), data=LongLS125.df, na.action = "na.exclude")
  tmp <- coef(summary(LMERlist_inter_Early[[i]]))
  tmpframe <- data.frame(inter_p = tmp[40], inter_Tstat = tmp[32], lausanne125label=labels_FWEmask[i])
  LMERslopes_inter_Early <- rbind(LMERslopes_inter_Early, tmpframe)
}

##P-value correcting##
LMERslopes_inter_Early <- LMERslopes_inter_Early[order(LMERslopes_inter_Early$inter_p),]
LMERslopes_inter_Early$FDRp <- p.adjust(LMERslopes_inter_Early$inter_p, 'fdr')
LMERslopes_inter_Early$FWEp <- p.adjust(LMERslopes_inter_Early$inter_p, 'bonferroni')

##Total LEQ##

##Initialize dataframes##
LMERlist_inter_Total <- vector(length=length(labels), mode='list')
LMERslopes_inter_Total <- data.frame()

##for loop generating model for each LS125 label##
for(i in 1:length(labels_FWEmask)) {
  LMERlist_inter_Total[[i]] <- lmer(as.formula(paste('lausanne125', labels[i], 'thickness_mean~(1|id)+TimeYR*Total.LEQ+AgeatMRI_base+DxDurYRS+SexFormatted+T1ProtocolFormatted', sep="_")), data=LongLS125.df, na.action = "na.exclude")
  tmp <- coef(summary(LMERlist_inter_Total[[i]]))
  tmpframe <- data.frame(inter_p = tmp[40], inter_Tstat = tmp[32], lausanne125label=labels_FWEmask[i])
  LMERslopes_inter_Total <- rbind(LMERslopes_inter_Total, tmpframe)
}

##P-value correcting##
LMERslopes_inter_Total <- LMERslopes_inter_Total[order(LMERslopes_inter_Total$inter_p),]
LMERslopes_inter_Total$FDRp <- p.adjust(LMERslopes_inter_Total$inter_p, 'fdr')
LMERslopes_inter_Total$FWEp <- p.adjust(LMERslopes_inter_Total$inter_p, 'bonferroni')

sig_LMER_inter_Total_unc005 <- subset(LMERslopes_inter_Total, LMERslopes_inter_Total$inter_p<0.005)
sig_LMER_inter_Total_unc005 <- data.frame(sig_LMER_inter_Total_unc005$lausanne125label, sig_LMER_inter_Total_unc005$inter_Tstat)
write.table(sig_LMER_inter_Total_unc005, file="Total_Interaction_SeqCovary_unc005.txt", row.names = FALSE, col.names = FALSE)

# Partial LEQ Regressions -------------------------------------------------

LateLEQmodel <- lm(Late.Life ~ Young.Adulthood + Mid.Life, data=CSLS125.df)
residuals <- LateLEQmodel$residuals
CSLS125_partial.df <- CSLS125.df
CSLS125_partial.df$Residuals <- residuals

LateLEQLMlist_partial <- vector(length=length(labels_FWEmask), mode='list')
LateLEQLMstats_partial <- data.frame()

##for loop that iterates thru the data frame of LS125 labels_FWEmask, LEQ scores, and covariates and creates a LM for LEQ predicting LS125 label (one for each LEQ stage)##
for(i in 1:length(labels_FWEmask)) {
  LateLEQLMlist_partial[[i]] <- lm(as.formula(paste('lausanne125', labels_FWEmask[i], 'thickness_mean~Residuals+AgeatMRI+SexFormatted+DxDurYRS+T1ProtocolFormatted', sep="_")), data=CSLS125_partial.df)
  tmp <- (LateLEQLMlist_partial[[i]])
  summtmp <- summary(tmp)
  coef <- coef(tmp)
  summcoef <- coef(summtmp)
  beta <- lm.beta(LateLEQLMlist_partial[[i]])$standardized.coefficients[2][[1]]
  adjrsq <- summtmp$adj.r.squared
  tmpframe <- data.frame(lausanne125label=labels_FWEmask[i], Intercept=summcoef[1], ResidualSlope=summcoef[2], ResidualTstat=summcoef[14], ResidualPval=summcoef[20],Beta=beta, Adj_Rsquared=adjrsq)
  LateLEQLMstats_partial <- rbind(LateLEQLMstats_partial, tmpframe)
}

##order by uncorrected p value##
LateLEQLMstats_partial <- LateLEQLMstats_partial[order(LateLEQLMstats_partial$ResidualPval),]

##FDR correction##
LateLEQLMstats_partial$FDRp<-p.adjust(LateLEQLMstats_partial$ResidualPval,'fdr')
LateLEQLMstats_partial$FWEp <- p.adjust(LateLEQLMstats_partial$ResidualPval, 'bonferroni')

##merge ls125 names into table##
LateLEQLMstats_partial <- merge(ls125_names.df, LateLEQLMstats_partial, by="lausanne125label")
LateLEQLMstats_partial <- LateLEQLMstats_partial[order(LateLEQLMstats_partial$ResidualPval),]

#output significant image##
sigLatestats_partial_unc005 <- subset(LateLEQLMstats_partial, LateLEQLMstats_partial$ResidualPval<0.005)

##new rendering format for Chris' Lausanne Render Script##
sigLatestats_partial_unc005 <- data.frame(sigLatestats_partial_unc005$lausanne125label, sigLatestats_partial_unc005$ResidualTstat)
write.table(sigLatestats_partial_unc005, file="CS_LateLEQ_Residual_unc005.txt", col.names = FALSE, row.names = FALSE)


##LONGITUDINAL -- SAME LME EXTRACTED SLOPES SO NO NEED TO RUN NEW LME##

##generate residuals##
LateLEQModel_long <- lm(Late.Life ~ Young.Adulthood + Mid.Life, data=SlopesforTesting)
residuals_long <- LateLEQModel_long$residuals
SlopesforTesting_partial <- SlopesforTesting
SlopesforTesting_partial$Residuals <- residuals_long

##initialize dataframes##
LateLongCorList_partial <- vector(length=length(labels_FWEmask), mode='list')
LateLongLMstats_partial <- data.frame()

##relate all LS125 region rate of change to late life LEQ##
for(i in 1:length(labels_FWEmask)){
  LateLongCorList_partial[[i]] <- lm(as.formula(paste('lausanne125', labels_FWEmask[i], 'slope~Residuals', sep="_")), data=SlopesforTesting_partial)
  tmp <- (LateLongCorList_partial[[i]])
  summcoef <- coef(summary(tmp))
  beta <- lm.beta(LateLongCorList_partial[[i]])$standardized.coefficients[2][[1]]
  adjrsq <- summary(tmp)$adj.r.squared
  tmpframe <- data.frame(lausanne125label=labels_FWEmask[i],ResidualP=summcoef[8], ResidualSlope=summcoef[2], ResidualTstat=summcoef[6],Beta=beta,Adj_Rsquared=adjrsq)
  LateLongLMstats_partial <- rbind(LateLongLMstats_partial, tmpframe)
}

##order dataframe and p-value correction##
LateLongLMstats_partial <- LateLongLMstats_partial[order(LateLongLMstats_partial$ResidualP),]
LateLongLMstats_partial$FDRp <- p.adjust(LateLongLMstats_partial$ResidualP, 'fdr')
LateLongLMstats_partial$FWEp <- p.adjust(LateLongLMstats_partial$ResidualP, 'bonferroni')

##merge ls125 names into table##
LateLongLMstats_partial <- merge(ls125_names.df, LateLongLMstats_partial, by="lausanne125label")
LateLongLMstats_partial <- LateLongLMstats_partial[order(LateLongLMstats_partial$ResidualP),]

unc005_LateLongLMstats_partial <- subset(LateLongLMstats_partial, LateLongLMstats_partial$ResidualP<0.005)
unc005_LateLongLMstats_partial <- data.frame(unc005_LateLongLMstats_partial$lausanne125label, unc005_LateLongLMstats_partial$ResidualTstat)
write.table(unc005_LateLongLMstats_partial, file="unc005_LateLEQ_ResidualLong.txt", col.names = FALSE, row.names = FALSE)  

