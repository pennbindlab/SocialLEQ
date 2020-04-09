##Generate IDTPs from INDD Query##

library(readxl)
library(MatchIt)

##read in data##
bvFTDCS.df <- read_excel("Updated_Cohort_3-17.xlsx", sheet="bvFTD-CS-Updated")
AllControls.df <- read_excel("Updated_Cohort_3-17.xlsx", sheet="PotentialControls")
bvFTDLong.df <- read_excel("Updated_Cohort_3-17.xlsx", sheet="bvFTD-Long")

##format data for propensity score matching##
bvMatch.df <- data.frame(bvFTDCS.df$id, bvFTDCS.df$date, bvFTDCS.df$AgeatMRI, bvFTDCS.df$SexFormatted, bvFTDCS.df$Education)
colnames(bvMatch.df) <- c("id","date","AgeatMRI","SexFormatted","Education")
bvMatch.df$ClinicalPhenotype1 <- "bvFTD"

CTRLmatch.df <- data.frame(AllControls.df$id, AllControls.df$date, AllControls.df$AgeatMRI, AllControls.df$SexFormatted, AllControls.df$Education, AllControls.df$ClinicalPhenotype1)
colnames(CTRLmatch.df) <- c("id","date","AgeatMRI","SexFormatted","Education","ClinicalPhenotype1")

FullMatch.df <- rbind(bvMatch.df, CTRLmatch.df)
FullMatch.df$IsPatient <- FullMatch.df$ClinicalPhenotype1 %in% 'bvFTD'

##Propensity score match and subset to make control cohort##
pm <- MatchIt::matchit(IsPatient ~ AgeatMRI + SexFormatted + Education, ratio=1, data=FullMatch.df, method='nearest')
FullMatch.df$Include <- pm$weights>0 

BVvsCTRL.df <- subset(FullMatch.df, FullMatch.df$Include=="TRUE")

##Save workplace for QuANTs script##

save.image("IDTPsupdated.Rdata", version = 2)
