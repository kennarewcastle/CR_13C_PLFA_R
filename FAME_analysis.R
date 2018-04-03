# 13C PLFA FAME Separation
# March 30, 2018
# KER
# Separating 13C signal in FAME compounds according to rhizsophere, substrate, and tree treatments

# Read in data
masterDat<-read.csv("13C_PLFA_master.csv")
#sampleNames<-data.frame(ID=seq(from=1,to=78),Sample_Name=unique(masterDat$Identifier.1))
#write.csv(sampleNames,file="sampleNames.csv")

# Obtain a list of FAME compounds identified in this analysis
# FAMEs<-unique(masterDat$Compounds.Result.Name)
# FAMEs # 60 different FAME compounds, 1 category for unknown FAME
# write.csv(FAMEs,file="FAME_compounds.csv")

# cleanDat<-masterDat[masterDat$Compounds.Result.Name!="",] # Removed compounds that are not biomarkers
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="10:00",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="11:00",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="13:00",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="2OH 12:0",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="3OH 13:0",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="14:01",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="3OH 14:0",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="17:01",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="2OH 16:0",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="18:3n6",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="trans-18:1n9",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="20:02",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="20:00",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="20:01",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="22:6n3",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="22:00",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="23:00",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="24:00:00",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="15:01",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="2OH 14:0",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="20:01",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="22:2 and 22:1n9",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="24:01:00",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="2OH 10:1",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="10:0",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="14:1",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="17:1",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="2O:0",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="22:0",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="24:0",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="20:1",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="21:00",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="11:0",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="13:0",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="20:2",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="21:0",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="23:0",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="15:1",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="24:1",]
# cleanDat<-cleanDat[cleanDat$Compounds.Result.Name!="20:0",]

#### Create a csv with only the meaningful compounds for each sample

#write.csv(cleanDat,file="CleanFAME_dat.csv")

#### Create function to calculate percent C in all fame compounds.

##########################################################################################################
# FUNCTION: percentC
# Function that calculates C % by mass for PLFA FAME compounds
# input: dat = dataframe of with molecular formula for each compound broken down into C,H,O columns
# output: perC = dataframe listing C % for each compound
#--------------------------------------------------------------------------------------------------------

percentC<-function(data=NULL) {
  if(is.null(data)) { # Minimalist code for default data frame
    comp_name<-seq(from=1,to=10)
    C<-sample(x=1:30,size=10,replace=TRUE)
    H<-sample(x=1:50,size=10,replace=TRUE)
    O<-rep(2,times=10)
    data<-data.frame(comp_name,C,H,O)
  }
  
  N<-length(data$C)
  perC<-c()
  
  for(i in 1:N) {
    C<-data$C[i] * 12.011
    H<-data$H[i] * 1.008
    O<-data$O[i] * 15.999
    molec_mass<-(C+H+O)
    perC[i]<-C/molec_mass
  }
  
  Cout<-data.frame(data,perC)
  return(Cout)
}

#_________________________________________________

#percentC()

fame_molec<-read.csv("fame_percent_C.csv")

fameCog<-percentC(data=fame_molec) # Data frame with percent_C calculations in last column

fameC<-data.frame(fameCog[,1:5],percent_C=fameCog[,7]) # Remove weird blank column

##### Creating dataframe that includes only variable relevant for analysis

cleanDat<-read.csv("CleanFAME_dat.csv") # This dataset has been edited using RegEx in BBedit to make compound labels consistent (ex. 14:00 and 14:0 = 14:0)
unique(cleanDat$Compounds.Result.Name)

####### Open CleanFAME_dat.csv in BBedit to make changes to sample names.