# 13C PLFA FAME Separation
# March 30, 2018
# KER
# Separating 13C signal in FAME compounds according to rhizsophere, substrate, and tree treatments

# Read in data
masterDat<-read.csv("13C_PLFA_master.csv")

# Obtain a list of FAME compounds identified in this analysis
FAMEs<-unique(masterDat$Compounds.Result.Name)
FAMEs # 60 different FAME compounds, 1 category for unknown FAME
write.csv(FAMEs,file="FAME_compounds.csv")

cleanDat<-masterDat[masterDat$Compounds.Result.Name!="",] # Removed compounds that are not biomarkers
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="10:00",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="11:00",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="13:00",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="2OH 12:0",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="3OH 13:0",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="14:01",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="3OH 14:0",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="17:01",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="2OH 16:0",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="18:3n6",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="trans-18:1n9",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="20:02",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="20:00",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="20:01",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="22:6n3",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="22:00",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="23:00",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="24:00:00",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="15:01",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="2OH 14:0",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="20:01",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="22:2 and 22:1n9",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="24:01:00",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="2OH 10:1",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="10:0",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="14:1",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="17:1",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="2O:0",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="22:0",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="24:0",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="20:1",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="21:00",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="11:0",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="13:0",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="20:2",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="21:0",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="23:0",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="15:1",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="24:1",]
cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="20:0",]




