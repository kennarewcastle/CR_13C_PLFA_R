### Calculation of carbon use efficiency from 13CO2 and 13C PLFA
### KER
### Created on: 16 November 2021


# Load packages -----------------------------------------------------------

library(dplyr)

# Function to convert raw PLFA data to ug C in FAME -----------------------

#########################################################################################################
# FUNCTION: ug_C_FAME
# Appends column with Peak Area:C standard ratio for each compound to data frame, uses this value to calcualte ugC of each fame compound, and appends to final data set.
# input: data = data frame where column 2= SampleName, column 6 = Peak_Area
#        stdRatio =  2-column data frame, column 1= sample name, column 2= standard PeakArea:C ratio for #        each sample.
#        molecular_mass = data frame with molecular mass for each FAME compound
# output: Data frame above with ugC for each FAME compound per sample is included as the last column in the data frame
#--------------------------------------------------------------------------------------------------------

ug_C_FAME<-function(data,stdRatio) {
  N<-length(data[,1])
  Cvec<-rep(NA,times=N)
  stdVec<-rep(NA,times=N)
  totalMassC<-rep(NA,times=N)
  
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
  
  ##### Divides ugC by injection volume, multiplies by 300 ul hexane to give total mass (ug) C in FAME     compound 
  for (i in 1:N) {
    totalMassC[i]<-(outDat2[i,12]/outDat2[i,3])*300
  } # close third for loop
  
  outDat3<-data.frame(outDat2,Total_C_ug=round(totalMassC,digits=2))
  

  finalDat<-data.frame(data,outDat3)
  return(finalDat)
  
} # close function body

#--------------------------------------------------------------------------------------------------------

## Read in data
data<-read.csv("MASTER_DAT_W_SOIL.csv") 
data$InjectionVol_ul<-as.numeric(data$InjectionVol_ul)
dat_stdRatio<-read.csv("PeakAreaC_ratio_per_sample.csv")
C_data<-ug_C_FAME(data=data,stdRatio=dat_stdRatio)

