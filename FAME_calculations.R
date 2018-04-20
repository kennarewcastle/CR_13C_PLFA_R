# Calculations for the concentrations of FAME compounds in each sample.
# 18 April 2018
# KER

library(dplyr)

# Read in data frame that was created with append functions in FAME_data_prep.R
data<-read.csv("MASTER_DAT_W_SOIL.csv") 
data<-data[,2:11] # Get rid of extra numerical ID column
data$InjectionVol_ul<-as.numeric(data$InjectionVol_ul)
# Filter out the FAME internal standards (12:0 and 19:0) from each sample
std_dat<-filter(data,FAME_compound=="12:0" | FAME_compound=="19:0")

##########################################################################################################
# FUNCTION: ugC
# Calculates ug C present in FAME standard samples that were added at 25 ug/ul
# input: data = data frame with the following columns for STANDARD FAMES in this order: ID, SampleName, InjectionVol_ul, FAME_compound, Sample_Dilution, Peak_Area, d13, Atom_per_13C, FAME_per_C, Dry_Soil_Mass_g
# output: Dataframe with column containing ugC per sample in the last column
#---------------------------------------------------------------------------------------------------------

ugC<-function(data) {
  N<-length(data[,1])
  Cvec<-rep("NA",times=N)
  
  for (i in 1:N) {
      Cvec[i]<-(data[i,3]*25)*data[i,9]
  }
datOUT<-data.frame(data,ug_C=as.numeric(Cvec))
return(datOUT)
}

#---------------------------------------------------------------------------------------------------------

std_perc<-ugC(data=std_dat)
std_perc$Peak_Area<-as.numeric(std_perc$Peak_Area)

##########################################################################################################
# FUNCTION: PeakC
# Calculates the ratio of peak area to ug C for each internal standard FAME compound.
# input: data = data frame to which peak area:C ratio will be appended; column 6 = Peak Area, column 11 = ugC in sample
# output: data frame with peak area:C for each FAME compound in the last column
#---------------------------------------------------------------------------------------------------------

PeakC<-function(data) {
  N<-length(data[,1])
  ratioVec<-rep(NA,times=N)
  
  for (i in 1:N) {
    ratioVec[i]<-data[i,6]/data[i,11]
  }
  ratioVec<-round(ratioVec,digits=4)
  outDat<-data.frame(data,Peak_Area_C_ratio=ratioVec)
  return(outDat)
}

#---------------------------------------------------------------------------------------------------------

data_PeakC<-PeakC(data=std_perc)
write.csv(data_PeakC,file="Standards_w_PeakArea_C_ratio.csv")

##########################################################################################################
# FUNCTION: stdRatio
# Calculates the mean Peak Area:C ratio for each sample from the two internal standards (19:0 and 12:0)
# input: data = data frame where column 2 = SampleName, column 12 = Peak_Area_C_ratio
# output: 
#---------------------------------------------------------------------------------------------------------

function<-function() {
  return("testing...........FunctionName")
}
