library(ANTsR)
library(lme4)
library(reshape2)
library(lmerTest)
library(MatchIt)

##loading data, format data and generate label numbers from QuANTs output##
load("LS125thick.Rdata")


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

sum(CSLS125.df$SexFormatted)
sum(CTRL.df$SexFormatted)
summary(CSLS125.df$AgeatMRI)
summary(CTRL.df$AgeatMRI)
summary(CSLS125.df$DxDurYRS)
summary(CSLS125.df$Education)
summary(CTRL.df$Education)

ctrlmmse.df <- merge(CTRL.df, AllControls.df, by="id")
summary(ctrlmmse.df$MMSETotal)

##Wilcox Testing between bv and CTRL predictor variables##
wilcox.test(AgeatMRI~Group, data=BVvsCTRLLS125.df)
wilcox.test(SexFormatted~Group, data=BVvsCTRLLS125.df)
wilcox.test(Education~Group, data=BVvsCTRLLS125.df)

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

##FDR correction and re-order##
Diffstats$GroupP.FDR<-p.adjust(Diffstats$GroupP,'fdr')
Diffstats$FWEP <- p.adjust(Diffstats$GroupP, 'bonferroni')

##take only significant p values for output##
sigDiffstats_FWE <- subset(Diffstats, Diffstats$FWEP < 0.05)
sigDiffstats_FDR <- subset(Diffstats, Diffstats$GroupP.FDR < 0.05)

#new visualization method Chris's script##
sigDiffstats_FWE <- data.frame(sigDiffstats_FWE$lausanne125label, sigDiffstats_FWE$FWEP)
write.table(sigDiffstats_FWE, file="BVvsCTRLstats_FWE.txt", col.names = FALSE, row.names = FALSE)

sigDiffstats_FDR <- data.frame(sigDiffstats_FDR$lausanne125label, sigDiffstats_FDR$GroupP.FDR)
write.table(sigDiffstats_FDR, file="BVvsCTRLstats_FDR.txt", col.names = FALSE, row.names = FALSE)

##Early LEQ##

##empty list to be filled with results of lm. empty dataframe for stats##
EarlyLEQLMlist <- vector(length=length(labels), mode='list')
EarlyLEQLMstats <- data.frame()

##for loop that iterates thru the data frame of LS125 labels, LEQ scores, and covariates and creates a LM for LEQ predicting LS125 label (one for each LEQ stage)##
for(i in 1:length(labels)) {
  EarlyLEQLMlist[[i]] <- lm(as.formula(paste('lausanne125', labels[i], 'thickness_mean~Young.Adulthood+AgeatMRI+SexFormatted+DxDurYRS', sep="_")), data=CSLS125.df)
  tmp <- (EarlyLEQLMlist[[i]])
  summtmp <- summary(tmp)
  coef <- coef(tmp)
  summcoef <- coef(summtmp)
  tmpframe <- data.frame(intercept=coef[1], EarlyLEQSlope=coef[2], EarlyLEQTstat=summcoef[12], EarlyLEQP=summcoef[17], AgeSlope=coef[3], AgeTstat=summcoef[13], AgeP=summcoef[18], SexSlope=coef[4], SexTstat=summcoef[14], SexP=summcoef[19], DxDurSlope=coef[5], DxDurTstat=summcoef[15], DxDurP=summcoef[20], lausanne125label=labels[i])
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
MidLEQLMlist <- vector(length=length(labels), mode='list')
MidLEQLMstats <- data.frame()

##for loop that iterates thru the data frame of LS125 labels, LEQ scores, and covariates and creates a LM for LEQ predicting LS125 label (one for each LEQ stage)##
for(i in 1:length(labels)) {
  MidLEQLMlist[[i]] <- lm(as.formula(paste('lausanne125', labels[i], 'thickness_mean~Mid.Life+AgeatMRI+SexFormatted+DxDurYRS', sep="_")), data=CSLS125.df)
  tmp <- (MidLEQLMlist[[i]])
  summtmp <- summary(tmp)
  coef <- coef(tmp)
  summcoef <- coef(summtmp)
  tmpframe <- data.frame(intercept=coef[1], MidLEQSlope=coef[2], MidLEQTstat=summcoef[12], MidLEQP=summcoef[17], AgeSlope=coef[3], AgeTstat=summcoef[13], AgeP=summcoef[18], SexSlope=coef[4], SexTstat=summcoef[14], SexP=summcoef[19], DxDurSlope=coef[5], DxDurTstat=summcoef[15], DxDurP=summcoef[20], lausanne125label=labels[i])
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

##Late LEQ##

##empty list to be filled with results of lm. empty dataframe for stats##
LateLEQLMlist <- vector(length=length(labels), mode='list')
LateLEQLMstats <- data.frame()

##for loop that iterates thru the data frame of LS125 labels, LEQ scores, and covariates and creates a LM for LEQ predicting LS125 label (one for each LEQ stage)##
for(i in 1:length(labels)) {
  LateLEQLMlist[[i]] <- lm(as.formula(paste('lausanne125', labels[i], 'thickness_mean~Late.Life+AgeatMRI+SexFormatted+DxDurYRS', sep="_")), data=CSLS125.df)
  tmp <- (LateLEQLMlist[[i]])
  summtmp <- summary(tmp)
  coef <- coef(tmp)
  summcoef <- coef(summtmp)
  tmpframe <- data.frame(intercept=coef[1], LateLEQSlope=coef[2], LateLEQTstat=summcoef[12], LateLEQP=summcoef[17], AgeSlope=coef[3], AgeTstat=summcoef[13], AgeP=summcoef[18], SexSlope=coef[4], SexTstat=summcoef[14], SexP=summcoef[19], DxDurSlope=coef[5], DxDurTstat=summcoef[15], DxDurP=summcoef[20], lausanne125label=labels[i])
  LateLEQLMstats <- rbind(LateLEQLMstats, tmpframe)
}

##order by uncorrected p value##
LateLEQLMstats <- LateLEQLMstats[order(LateLEQLMstats$LateLEQP),]

##FDR correction##
LateLEQLMstats$FDRp<-p.adjust(LateLEQLMstats$LateLEQP,'fdr')
LateLEQLMstats$FWep <- p.adjust(LateLEQLMstats$LateLEQP, 'bonferroni')

#output significant image##
sigLatestats <- subset(LateLEQLMstats, LateLEQLMstats$FDRp<0.05)
sigLatestats <- subset(sigLatestats, sigLatestats$LateLEQSlope>0)

##new rendering format for Chris' Lausanne Render Script##
sigLatestats <- data.frame(sigLatestats$lausanne125label, sigLatestats$FDRp)
write.table(sigLatestats, file="LateCSStats.txt", col.names = FALSE, row.names = FALSE)

##Total LEQ##

TotalLEQLMlist <- vector(length=length(labels), mode='list')
TotalLeqLMstats <- data.frame()

##for loop that iterates thru the data frame of LS125 labels, LEQ scores, and covariates and creates a LM for LEQ predicting LS125 label (one for each LEQ stage)##
for(i in 1:length(labels)) {
  TotalLEQLMlist[[i]] <- lm(as.formula(paste('lausanne125', labels[i], 'thickness_mean~Total.LEQ+AgeatMRI+SexFormatted+DxDurYRS', sep="_")), data=CSLS125.df)
  tmp <- (TotalLEQLMlist[[i]])
  summtmp <- summary(tmp)
  coef <- coef(tmp)
  summcoef <- coef(summtmp)
  tmpframe <- data.frame(intercept=coef[1], TotalLEQSlope=coef[2], TotalLEQTstat=summcoef[12], TotalLEQP=summcoef[17], AgeSlope=coef[3], AgeTstat=summcoef[13], AgeP=summcoef[18], SexSlope=coef[4], SexTstat=summcoef[14], SexP=summcoef[19], DxDurSlope=coef[5], DxDurTstat=summcoef[15], DxDurP=summcoef[20], lausanne125label=labels[i])
  TotalLeqLMstats <- rbind(TotalLeqLMstats, tmpframe)
}

##order by uncorrected p value##
TotalLeqLMstats <- TotalLeqLMstats[order(TotalLeqLMstats$TotalLEQP),]

##FDR correction##
TotalLeqLMstats$FDRp<-p.adjust(TotalLeqLMstats$TotalLEQP,'fdr')
TotalLeqLMstats$FWep <- p.adjust(TotalLeqLMstats$TotalLEQP, 'bonferroni')

#output significant image##
sigTotalstats <- subset(TotalLeqLMstats, TotalLeqLMstats$FDRp<0.05)
sigTotalstats <- subset(sigTotalstats, sigTotalstats$TotalLEQSlope>0)

##new rendering format for Chris' Lausanne Render Script##
sigTotalstats <- data.frame(sigTotalstats$lausanne125label, sigTotalstats$FDRp)
write.table(sigTotalstats, file="TotalCSStats.txt", col.names = FALSE, row.names = FALSE)

##linear mixed effects model to determine how LEQ affects rate of decline##

##generate model and extract slope to relate rate of decline to LEQ##
library(lme4)

##empty list to be filled with results of lmer. empty dataframe for slopes##
LMERlist <- vector(length=length(labels), mode='list')
LMERslopes <- data.frame()
singtest <- data.frame()

##for loop generating model for each LS125 label##
for(i in 1:length(labels)) {
  LMERlist[[i]] <- lmer(as.formula(paste('lausanne125', labels[i], 'thickness_mean~BaseDxDurYR + (TimeYR|id) + SexFormatted + AgeatMRI', sep="_")), data=LongLS125.df, na.action = "na.exclude", control=lmerControl(check.nobs.vs.nRE="ignore"))
  sing <- isSingular(LMERlist[[i]])
  tempsing <- data.frame(label=paste('lausanne125',labels[i],sep="_"), sing=sing)
  singtest <- rbind(singtest,tempsing)
  tmp <- (LMERlist[[i]])
  placeholder <- coef(tmp)
  temp <- placeholder$id
  slopes <- setNames(cbind(rownames(temp), temp,row.names=NULL),c("id","Slopes","Intercept","NA1","NA2","NA3"))
  tmpframe <- data.frame(id=slopes$id, slope=slopes$Slopes, lausanne125label=paste('lausanne125',labels[i],'slope',sep="_"))
  LMERslopes <- rbind(LMERslopes, tmpframe)
}

##reformat for further testing and add demographic values from INDD##
SlopesforTesting <- dcast(LMERslopes, id~lausanne125label, value.var="slope")
SlopesforTesting <- merge(SlopesforTesting, tp1, by="id")

##relate slopes to LEQ##

##Early##

##empty list to fill with results of lm##
EarlyLongCorList <- vector(length=length(labels), mode='list')
EarlyLongLMstats <- data.frame()

##reEarly all LS125 region rate of change to Early life LEQ##
for(i in 1:length(labels)){
  EarlyLongCorList[[i]] <- lm(as.formula(paste('lausanne125', labels[i], 'slope~Young.Adulthood', sep="_")), data=SlopesforTesting)
  tmp <- (EarlyLongCorList[[i]])
  summcoef <- coef(summary(tmp))
  tmpframe <- data.frame(EarlyLEQP=summcoef[8], EarlyLEQSlope=summcoef[2], lausanne125label=labels[i])
  EarlyLongLMstats <- rbind(EarlyLongLMstats, tmpframe)
}

##order dataframe and p-value correction##
EarlyLongLMstats <- EarlyLongLMstats[order(EarlyLongLMstats$EarlyLEQP),]
EarlyLongLMstats$FDRp <- p.adjust(EarlyLongLMstats$EarlyLEQP, 'fdr')
EarlyLongLMstats$FWEp <- p.adjust(EarlyLongLMstats$EarlyLEQP, 'bonferroni')

##get significant regions##
sigEarlyLongLMstats <- subset(EarlyLongLMstats, EarlyLongLMstats$FDRp<0.05)
sigEarlyLongLMstats <- subset(sigEarlyLongLMstats, sigEarlyLongLMstats$EarlyLEQSlope>0)

##new rendering format for Chris' Lausanne Render Script##
sigEarlyLongLMstats <- data.frame(sigEarlyLongLMstats$lausanne125label, sigEarlyLongLMstats$FDRp)
write.table(sigEarlyLongLMstats, file="EarlyLongLMStats.txt", col.names = FALSE, row.names = FALSE)


##Mid##

##empty list to fill with results of lm##
MidLongCorList <- vector(length=length(labels), mode='list')
MidLongLMstats <- data.frame()

##reMid all LS125 region rate of change to Mid life LEQ##
for(i in 1:length(labels)){
  MidLongCorList[[i]] <- lm(as.formula(paste('lausanne125', labels[i], 'slope~Mid.Life', sep="_")), data=SlopesforTesting)
  tmp <- (MidLongCorList[[i]])
  summcoef <- coef(summary(tmp))
  tmpframe <- data.frame(MidLEQP=summcoef[8], MidLEQSlope=summcoef[2], lausanne125label=labels[i])
  MidLongLMstats <- rbind(MidLongLMstats, tmpframe)
}

##order dataframe and p-value correction##
MidLongLMstats <- MidLongLMstats[order(MidLongLMstats$MidLEQP),]
MidLongLMstats$FDRp <- p.adjust(MidLongLMstats$MidLEQP, 'fdr')
MidLongLMstats$FWEp <- p.adjust(MidLongLMstats$MidLEQP, 'bonferroni')

##get significant regions##
sigMidLongLMstats <- subset(MidLongLMstats, MidLongLMstats$FDRp<0.05)
sigMidLongLMstats <- subset(sigMidLongLMstats, sigMidLongLMstats$MidLEQSlope>0)

##new rendering format for Chris' Lausanne Render Script##
sigMidLongLMstats <- data.frame(sigMidLongLMstats$lausanne125label, sigMidLongLMstats$FDRp)
write.table(sigMidLongLMstats, file="MidLongLMStats.txt", col.names = FALSE, row.names = FALSE)


##Late##

##empty list to fill with results of lm##
LateLongCorList <- vector(length=length(labels), mode='list')
LateLongLMstats <- data.frame()

##relate all LS125 region rate of change to late life LEQ##
for(i in 1:length(labels)){
  LateLongCorList[[i]] <- lm(as.formula(paste('lausanne125', labels[i], 'slope~Late.Life', sep="_")), data=SlopesforTesting)
  tmp <- (LateLongCorList[[i]])
  summcoef <- coef(summary(tmp))
  tmpframe <- data.frame(LateLEQP=summcoef[8], LateLEQSlope=summcoef[2], lausanne125label=labels[i])
  LateLongLMstats <- rbind(LateLongLMstats, tmpframe)
}

##order dataframe and p-value correction##
LateLongLMstats <- LateLongLMstats[order(LateLongLMstats$LateLEQP),]
LateLongLMstats$FDRp <- p.adjust(LateLongLMstats$LateLEQP, 'fdr')
LateLongLMstats$FWEp <- p.adjust(LateLongLMstats$LateLEQP, 'bonferroni')

##get significant regions##
sigLateLongLMstats <- subset(LateLongLMstats, LateLongLMstats$FDRp<0.05)
sigLateLongLMstats <- subset(sigLateLongLMstats, sigLateLongLMstats$LateLEQSlope>0)

##new rendering format for Chris' Lausanne Render Script##
sigLateLongLMstats <- data.frame(sigLateLongLMstats$lausanne125label, sigLateLongLMstats$FDRp)
write.table(sigLateLongLMstats, file="LateLongLMStats.txt", col.names = FALSE, row.names = FALSE)


##Total##

##empty list to fill with results of lm##
TotalLongCorList <- vector(length=length(labels), mode='list')
TotalLongLMstats <- data.frame()

##reTotal all LS125 region rate of change to Total life LEQ##
for(i in 1:length(labels)){
  TotalLongCorList[[i]] <- lm(as.formula(paste('lausanne125', labels[i], 'slope~Total.LEQ', sep="_")), data=SlopesforTesting)
  tmp <- (TotalLongCorList[[i]])
  summcoef <- coef(summary(tmp))
  tmpframe <- data.frame(TotalLEQP=summcoef[8], TotalLEQSlope=summcoef[2], lausanne125label=labels[i])
  TotalLongLMstats <- rbind(TotalLongLMstats, tmpframe)
}

##order dataframe and p-value correction##
TotalLongLMstats <- TotalLongLMstats[order(TotalLongLMstats$TotalLEQP),]
TotalLongLMstats$FDRp <- p.adjust(TotalLongLMstats$TotalLEQP, 'fdr')
TotalLongLMstats$FWEp <- p.adjust(TotalLongLMstats$TotalLEQP, 'bonferroni')

##get significant regions##
sigTotalLongLMstats <- subset(TotalLongLMstats, TotalLongLMstats$FDRp<0.05)
sigTotalLongLMstats <- subset(sigTotalLongLMstats, sigTotalLongLMstats$TotalLEQSlope>0)

##new rendering format for Chris' Lausanne Render Script##
sigTotalLongLMstats <- data.frame(sigTotalLongLMstats$lausanne125label, sigTotalLongLMstats$FDRp)
write.table(sigTotalLongLMstats, file="TotalLongLMStats.txt", col.names = FALSE, row.names = FALSE)




