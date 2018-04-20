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
write.csv(data_PeakC,file="Standards_w_PeakArea_C_ratio.csv",row.names=FALSE) # row.names=FALSE keeps extra row index from being added

##########################################################################################################
# FUNCTION: stdRatio
# Calculates the mean Peak Area:C ratio for each sample from the two internal standards (19:0 and 12:0)
# input: data = data frame where column 2 = SampleName, column 12 = Peak_Area_C_ratio
# output: New data frame with an average standard PeakArea:C ratio for each sample
#---------------------------------------------------------------------------------------------------------
library(dplyr)

stdRatio<-function(data) {
  N<-length(data[,1])
  ratioVec<-rep(NA,times=N)
  
  for (i in 1:N) {
   if (i == 1) {
     x<-c(data[1,12],data[2,12])
     ratioVec[i]<-mean(x)
      } else {
    if (data[i,2] == data[i-1,2]) {
      ratioVec[i]<-ratioVec[i-1]
      }
    if (data[i,2] != data[i-1,2]) {
      x<-c(data[i,12],data[i+1,12])
      ratioVec[i]<-mean(x)
      } # close final if statement
    } # close else statement
  } # close for loop
    
  pre_outDat<-data.frame(SampleName=data[,2],Mean_PeakArea_C_Ratio=as.numeric(ratioVec))
  outDat<-pre_outDat[!duplicated(pre_outDat[,'SampleName']),] # Pulls out 1 average standard ratio for each sample.
  return(outDat)
  
} # close function body

#--------------------------------------------------------------------------------------------------------

data_PeakC<-read.csv("Standards_w_PeakArea_C_ratio.csv")

dat_stdRatio<-stdRatio(data=data_PeakC)
write.csv(dat_stdRatio,file="PeakAreaC_ratio_per_sample.csv",row.names=FALSE)

#########################################################################################################
# FUNCTION: nmolFAME
# Appends column with Peak Area:C standard ratio for each compound to data frame, uses this value to calcualte ugC of each fame compound, and appends to final data set.
# input: data = data frame where column 2= SampleName, column 6 = Peak_Area
#        stdRatio =  2-column data frame, column 1= sample name, column 2= standard PeakArea:C ratio for #        each sample.
# output: Data frame above with ugC for each FAME compound per sample is included as the last column in the data frame
#--------------------------------------------------------------------------------------------------------

nmolFAME<-function(data,stdRatio) {
  N<-length(data[,1])
  Cvec<-rep(NA,times=N)
  stdVec<-rep(NA,times=N)
  totalMassC<-rep(NA,times=N)
  ugFAME<-rep(NA,times=N)
  
  for (i in 1:N) {
    if(data[i,2]=="01.05.1"){
      stdVec[i]<-stdRatio[1,2]
    }
    if(data[i,2]=="01.13.1"){
      stdVec[i]<-stdRatio[2,2]
    }
    if(data[i,2]=="01.13.2"){
      stdVec[i]<-stdRatio[3,2]
    }
    if(data[i,2]=="01.13.3"){
      stdVec[i]<-stdRatio[4,2]
    }
    if(data[i,2]=="01.05.3"){
      stdVec[i]<-stdRatio[5,2]
    }
    if(data[i,2]=="10.02.2"){
      stdVec[i]<-stdRatio[6,2]
    }
    if(data[i,2]=="10.02.1"){
      stdVec[i]<-stdRatio[7,2]
    }
    if(data[i,2]=="10.02.3"){
      stdVec[i]<-stdRatio[8,2]
    }
    if(data[i,2]=="11.06.2"){
      stdVec[i]<-stdRatio[9,2]
    }
    if(data[i,2]=="11.10.1"){
      stdVec[i]<-stdRatio[10,2]
    }
    if(data[i,2]=="11.10.3"){
      stdVec[i]<-stdRatio[11,2]
    }
    if(data[i,2]=="11.06.3"){
      stdVec[i]<-stdRatio[12,2]
    }
    if(data[i,2]=="12.09.1"){
      stdVec[i]<-stdRatio[13,2]
    }
    if(data[i,2]=="12.09.2"){
      stdVec[i]<-stdRatio[14,2]
    }
    if(data[i,2]=="12.01.2"){
      stdVec[i]<-stdRatio[15,2]
    }
    if(data[i,2]=="12.01.3"){
      stdVec[i]<-stdRatio[16,2]
    }
    if(data[i,2]=="12.05.2"){
      stdVec[i]<-stdRatio[17,2]
    }
    if(data[i,2]=="12.05.3"){
      stdVec[i]<-stdRatio[18,2]
    }
    if(data[i,2]=="12.09.3"){
      stdVec[i]<-stdRatio[19,2]
    }
    if(data[i,2]=="14.10.3"){
      stdVec[i]<-stdRatio[20,2]
    }
    if(data[i,2]=="14.10.1"){
      stdVec[i]<-stdRatio[21,2]
    }
    if(data[i,2]=="14.10.2"){
      stdVec[i]<-stdRatio[22,2]
    }
    if(data[i,2]=="15.16.2"){
      stdVec[i]<-stdRatio[23,2]
    }
    if(data[i,2]=="17.01.1"){
      stdVec[i]<-stdRatio[24,2]
    }
    if(data[i,2]=="17.01.3"){
      stdVec[i]<-stdRatio[25,2]
    }
    if(data[i,2]=="17.06.1"){
      stdVec[i]<-stdRatio[26,2]
    }
    if(data[i,2]=="17.06.2"){
      stdVec[i]<-stdRatio[27,2]
    }
    if(data[i,2]=="17.01.1b"){
      stdVec[i]<-stdRatio[28,2]
    }
    if(data[i,2]=="02.08.1"){
      stdVec[i]<-stdRatio[35,2]
    }
    if(data[i,2]=="02.08.2"){
      stdVec[i]<-stdRatio[36,2]
    }
    if(data[i,2]=="02.08.3"){
      stdVec[i]<-stdRatio[37,2]
    }
    if(data[i,2]=="02.10.1"){
      stdVec[i]<-stdRatio[33,2]
    }
    if(data[i,2]=="02.10.3"){
      stdVec[i]<-stdRatio[34,2]
    }
    if(data[i,2]=="03.07.2"){
      stdVec[i]<-stdRatio[47,2]
    }
    if(data[i,2]=="03.07.3"){
      stdVec[i]<-stdRatio[48,2]
    }
    if(data[i,2]=="05.03.1"){
      stdVec[i]<-stdRatio[61,2]
    }
    if(data[i,2]=="05.03.2"){
      stdVec[i]<-stdRatio[62,2]
    }
    if(data[i,2]=="05.03.3"){
      stdVec[i]<-stdRatio[63,2]
    }
    if(data[i,2]=="05.09.2"){
      stdVec[i]<-stdRatio[64,2]
    }
    if(data[i,2]=="05.09.3"){
      stdVec[i]<-stdRatio[65,2]
    }
    if(data[i,2]=="06.06.1"){
      stdVec[i]<-stdRatio[74,2]
    }
    if(data[i,2]=="06.06.2"){
      stdVec[i]<-stdRatio[75,2]
    }
    if(data[i,2]=="06.06.3"){
      stdVec[i]<-stdRatio[76,2]
    }
    if(data[i,2]=="17.06.3"){
      stdVec[i]<-stdRatio[29,2]
    }
    if(data[i,2]=="18.05.1"){
      stdVec[i]<-stdRatio[30,2]
    }
    if(data[i,2]=="18.05.2"){
      stdVec[i]<-stdRatio[31,2]
    }
    if(data[i,2]=="18.05.3"){
      stdVec[i]<-stdRatio[32,2]
    }
    if(data[i,2]=="20.04.1"){
      stdVec[i]<-stdRatio[38,2]
    }
    if(data[i,2]=="20.04.2"){
      stdVec[i]<-stdRatio[39,2]
    }
    if(data[i,2]=="22.07.2"){
      stdVec[i]<-stdRatio[40,2]
    }
    if(data[i,2]=="22.07.3"){
      stdVec[i]<-stdRatio[41,2]
    }
    if(data[i,2]=="24.02.1"){
      stdVec[i]<-stdRatio[42,2]
    }
    if(data[i,2]=="24.02.2"){
      stdVec[i]<-stdRatio[43,2]
    }
    if(data[i,2]=="24.02.3"){
      stdVec[i]<-stdRatio[44,2]
    }
    if(data[i,2]=="28.01.2"){
      stdVec[i]<-stdRatio[45,2]
    }
    if(data[i,2]=="28.01.3"){
      stdVec[i]<-stdRatio[46,2]
    }
    if(data[i,2]=="34.06.1"){
      stdVec[i]<-stdRatio[49,2]
    }
    if(data[i,2]=="34.06.3"){
      stdVec[i]<-stdRatio[50,2]
    }
    if(data[i,2]=="37.04.1"){
      stdVec[i]<-stdRatio[51,2]
    }
    if(data[i,2]=="37.04.2"){
      stdVec[i]<-stdRatio[52,2]
    }
    if(data[i,2]=="37.04.3"){
      stdVec[i]<-stdRatio[53,2]
    }
    if(data[i,2]=="43.02.1"){
      stdVec[i]<-stdRatio[54,2]
    }
    if(data[i,2]=="43.02.2"){
      stdVec[i]<-stdRatio[55,2]
    }
    if(data[i,2]=="43.02.3"){
      stdVec[i]<-stdRatio[56,2]
    }
    if(data[i,2]=="46.03.2"){
      stdVec[i]<-stdRatio[57,2]
    }
    if(data[i,2]=="46.03.3"){
      stdVec[i]<-stdRatio[58,2]
    }
    if(data[i,2]=="48.01.1"){
      stdVec[i]<-stdRatio[59,2]
    }
    if(data[i,2]=="48.01.3"){
      stdVec[i]<-stdRatio[60,2]
    }
    if(data[i,2]=="52.04.1"){
      stdVec[i]<-stdRatio[67,2]
    }
    if(data[i,2]=="52.04.2"){
      stdVec[i]<-stdRatio[66,2]
    }
    if(data[i,2]=="52.04.3"){
      stdVec[i]<-stdRatio[68,2]
    }
    if(data[i,2]=="55.11.1"){
      stdVec[i]<-stdRatio[69,2]
    }
    if(data[i,2]=="55.11.2"){
      stdVec[i]<-stdRatio[70,2]
    }
    if(data[i,2]=="58.05.1"){
      stdVec[i]<-stdRatio[73,2]
    }
    if(data[i,2]=="58.05.2"){
      stdVec[i]<-stdRatio[71,2]
    }
    if(data[i,2]=="58.05.3"){
      stdVec[i]<-stdRatio[72,2]
    }
    if(data[i,2]=="63.04.2"){
      stdVec[i]<-stdRatio[77,2]
    } # close final if statement
  } # close first for loop
  
  outDat1<-data.frame(data,StandardRatio=stdVec)

###### Calculates ugC present in each compound from Peak Area
  for (i in 1:N) {
    Cvec[i]<-outDat1[i,6]/outDat1[i,11]
    Cvec[i]<-round(Cvec[i],digits=4)
  } # close second for loop
  
  outDat2<-data.frame(outDat1,ugC_FAME_compound=Cvec)

##### Divides ugC by injection volume, multiplies by 300 ul hexane to give total mass (ug) C in FAME compound 
  for (i in 1:N) {
    totalMassC[i]<-(outDat2[i,12]/outDat2[i,3])*300
  } # close third for loop
    
  outDat3<-data.frame(outDat2,Total_C_ug=round(totalMassC,digits=2))
  
##### Divides by the percent carbon in each FAME compound to yield ug FAME compound in total sample
  for (i in 1:N) {
    ugFAME[i]<-outDat3[i,13]/outDat3[i,9]
  } # close fourth for loop
  
  outDat4<-data.frame(outDat3,Total_Mass_FAME_ug=round(ugFAME,digits=2))
  return(outDat4)
  
} # close function body

nmolFAME(data=data,stdRatio=dat_stdRatio)

#### Data to run the nmolFAAME function
data<-read.csv("MASTER_DAT_W_SOIL.csv") 
data<-data[,2:11] # Get rid of extra numerical ID column
data$InjectionVol_ul<-as.numeric(data$InjectionVol_ul)

dat_stdRatio<-read.csv("PeakAreaC_ratio_per_sample.csv")

               