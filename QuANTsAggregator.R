library(QuANTs)

load("IDTPsupdated.Rdata")

path="/data/grossman/pipedream2018//crossSectional/antsct/"

CSIDs <- as.character(bvFTDCS.df$id)
CSTPs <- as.character(bvFTDCS.df$date)

LongIDs <- as.character(bvFTDLong.df$id)
LongTPs <- as.character(bvFTDLong.df$date)

BVvCTRLIDs <- as.character(BVvsCTRL.df$id)
BVvCTRLTPs <- as.character(BVvsCTRL.df$date)

CSLS125.df <- getQuants(path=path, id=CSIDs, date=CSTPs, system="lausanne125", measure="thickness", metric="mean", as.wide=T)
LongLS125.df <- getQuants(path=path, id=LongIDs, date=LongTPs, system="lausanne125", measure="thickness", metric="mean", as.wide=T)
BVvsCTRLLS125.df <- getQuants(path=path, id=BVvCTRLIDs, date=BVvCTRLTPs, system="lausanne125", measure="thickness", metric="mean", as.wide=T)

save.image("LS125thick.Rdata")
