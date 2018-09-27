# Calculations for the concentrations of FAME compounds in each sample.
# 18 April 2018
# Last modified: 
# KER

library(dplyr)

# Read in data frame that was created with append functions in FAME_data_prep.R
data<-read.csv("MASTER_DAT_W_SOIL.csv") 
data$InjectionVol_ul<-as.numeric(data$InjectionVol_ul)

# Filter out the FAME internal standards (12:0 and 19:0) from each sample
std_dat<-filter(data,FAME_compound=="12:00" | FAME_compound=="19:00") # change this back to single 0s after editing compound names

#########################################################################################################
# FUNCTION: ugC
# Calculates ug C present in FAME standard samples that were added at 25 ug/ul
# input: data = data frame with the following columns for STANDARD FAMES in this order: ID, SampleName, InjectionVol_ul, FAME_compound, Sample_Dilution, Peak_Area, d13, Atom_per_13C, FAME_per_C, Dry_Soil_Mass_g
# output: Dataframe with column containing ugC per sample in the last column
#--------------------------------------------------------------------------------------------------------

ugC<-function(data) {
  N<-length(data[,1])
  Cvec<-rep("NA",times=N)
  
  for (i in 1:N) {
      Cvec[i]<-(data[i,3]*25)*data[i,9]
  }
datOUT<-data.frame(data,ug_C=as.numeric(Cvec))
return(datOUT)
}

#--------------------------------------------------------------------------------------------------------

std_perc<-ugC(data=std_dat)
std_perc$Peak_Area<-as.numeric(std_perc$Peak_Area)

#########################################################################################################
# FUNCTION: PeakC
# Calculates the ratio of peak area to ug C for each internal standard FAME compound.
# input: data = data frame to which peak area:C ratio will be appended; column 6 = Peak Area, column 11 = ugC in sample
# output: data frame with peak area:C for each FAME compound in the last column
#--------------------------------------------------------------------------------------------------------

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

#--------------------------------------------------------------------------------------------------------

data_PeakC<-PeakC(data=std_perc)
write.csv(data_PeakC,file="Standards_w_PeakArea_C_ratio.csv",row.names=FALSE) # row.names=FALSE keeps extra row index from being added

#########################################################################################################
# FUNCTION: stdRatio
# Calculates the mean Peak Area:C ratio for each sample from the two internal standards (19:0 and 12:0)
# input: data = data frame where column 2 = SampleName, column 12 = Peak_Area_C_ratio
# output: New data frame with an average standard PeakArea:C ratio for each sample
#--------------------------------------------------------------------------------------------------------
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
# FUNCTION: MolecMass
# Calculates molecular mass of each FAME compound, prints data frame
# input: data = data frame with the number of C, H, and O atoms in each compound in individual columns
# output: original data frame with percent_C and molec_mass calculated and appended as columns at the end
#--------------------------------------------------------------------------------------------------------

MolecMass<-function(data) {
  N<-length(data[,1])
  molec_mass<-rep(NA,times=N)
  perC<-rep(NA,times=N)
  
  for(i in 1:N) {
    C<-data$C[i] * 12.011
    H<-data$H[i] * 1.008
    O<-data$O[i] * 15.999
    molec_mass[i]<-(C+H+O)
    perC[i]<-C/molec_mass[i]
  }
  
  outDat<-data.frame(data,percent_C=perC,molecular_mass=molec_mass)
  return(outDat)
}

#--------------------------------------------------------------------------------------------------------

FAME_elements<-read.csv("fame_percent_C.csv")
FAME_for_csv<-MolecMass(data=FAME_elements)
write.csv(FAME_for_csv,file="FAME_molecular_mass_data.csv",row.names=FALSE)

#########################################################################################################
# FUNCTION: nmolFAME
# Appends column with Peak Area:C standard ratio for each compound to data frame, uses this value to calcualte ugC of each fame compound, and appends to final data set.
# input: data = data frame where column 2= SampleName, column 6 = Peak_Area
#        stdRatio =  2-column data frame, column 1= sample name, column 2= standard PeakArea:C ratio for #        each sample.
#        molecular_mass = data frame with molecular mass for each FAME compound
# output: Data frame above with ugC for each FAME compound per sample is included as the last column in the data frame
#--------------------------------------------------------------------------------------------------------

nmolFAME<-function(data,stdRatio,molecular_mass) {
  N<-length(data[,1])
  Cvec<-rep(NA,times=N)
  stdVec<-rep(NA,times=N)
  totalMassC<-rep(NA,times=N)
  ugFAME<-rep(NA,times=N)
  molec_mass_vec<-rep(NA,times=N)
  nmolVec<-rep(NA,times=N)
  FAMEsoil<-rep(NA,times=N)
  
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
  
  ##### Divides by the percent carbon in each FAME compound to yield ug FAME compound in total sample
  for (i in 1:N) {
    ugFAME[i]<-outDat3[i,13]/outDat3[i,9]
  } # close fourth for loop
  
  outDat4<-data.frame(outDat3,Total_Mass_FAME_ug=round(ugFAME,digits=2))

  ##### Appends molecular mass of each FAME compound to data frame.
  for (i in 1:N) {
    if(outDat4[i,4]=="12:0"){
      molec_mass_vec[i]<-round(molecular_mass[1,7],digits=4)
    }
    if(outDat4[i,4]=="14:0"){
      molec_mass_vec[i]<-round(molecular_mass[2,7],digits=4)
    }
    if(outDat4[i,4]=="a15:0"){
      molec_mass_vec[i]<-round(molecular_mass[3,7],digits=4)
    }
    if(outDat4[i,4]=="15:0"){
      molec_mass_vec[i]<-round(molecular_mass[4,7],digits=4)
    }
    if(outDat4[i,4]=="16:1n9/i16:0"){
      molec_mass_vec[i]<-round(molecular_mass[5,7],digits=4)
    }
    if(outDat4[i,4]=="16:0"){
      molec_mass_vec[i]<-round(molecular_mass[6,7],digits=4)
    }
    if(outDat4[i,4]=="i17:0"){
      molec_mass_vec[i]<-round(molecular_mass[7,7],digits=4)
    }
    if(outDat4[i,4]=="17:0cy"){
      molec_mass_vec[i]<-round(molecular_mass[8,7],digits=4)
    }
    if(outDat4[i,4]=="17:0"){
      molec_mass_vec[i]<-round(molecular_mass[9,7],digits=4)
    }
    if(outDat4[i,4]=="18:0"){
      molec_mass_vec[i]<-round(molecular_mass[10,7],digits=4)
    }
    if(outDat4[i,4]=="19:0cy"){
      molec_mass_vec[i]<-round(molecular_mass[11,7],digits=4)
    }
    if(outDat4[i,4]=="19:0"){
      molec_mass_vec[i]<-round(molecular_mass[12,7],digits=4)
    }
    if(outDat4[i,4]=="20:4n6"){
      molec_mass_vec[i]<-round(molecular_mass[13,7],digits=4)
    }
    if(outDat4[i,4]=="20:5n3"){
      molec_mass_vec[i]<-round(molecular_mass[14,7],digits=4)
    }
    if(outDat4[i,4]=="18:2n9,12 and cis18:1n9 and 18:3n3"){
      molec_mass_vec[i]<-round(molecular_mass[15,7],digits=4)
    }
    if(outDat4[i,4]=="18:2n9,12 and cis-18:1n9 and 183n3"){
      molec_mass_vec[i]<-round(molecular_mass[15,7],digits=4)
    }
    if(outDat4[i,4]=="trans-18:1n9"){
      molec_mass_vec[i]<-round(molecular_mass[16,7],digits=4)
    }
    if(outDat4[i,4]=="18:3n6"){
      molec_mass_vec[i]<-round(molecular_mass[17,7],digits=4)
    }
    if(outDat4[i,4]=="16:1w5"){
      molec_mass_vec[i]<-round(molecular_mass[18,7],digits=4)
    }
  } # close fifth for loop
  
  outDat5<-data.frame(outDat4,FAME_molecular_mass=molec_mass_vec)
  
  ##### Calculates nmol FAME compound per sample using molecular mass of FAME compound
  for (i in 1:N) {
    g<-outDat5[i,14]/1000000 # ug FAME to g FAME
    moles<-g/outDat5[i,15] # converts g to moles using molecular mass
    nmolVec[i]<-moles * 10^9 # converts moles to nanomoles
  } # close sixth for loop

  outDat6<-data.frame(outDat5,nmol_FAME=nmolVec)
 
  ##### Calculates nmol FAME per g dry soil in sample
  for (i in 1:N) {
    FAMEsoil[i]<-outDat6[i,16]/outDat6[i,10] # Divide nmol FAME by g dry soil in each sample
  } # close seventh for loop
    
  finalDat<-data.frame(data,nmol_FAME_per_g_soil=round(FAMEsoil,digits=2))
  return(finalDat)
  
} # close function body

#--------------------------------------------------------------------------------------------------------

#### Data to run the nmolFAME function
data<-read.csv("MASTER_DAT_W_SOIL.csv") 
data$InjectionVol_ul<-as.numeric(data$InjectionVol_ul)

dat_stdRatio<-read.csv("PeakAreaC_ratio_per_sample.csv")
dat_molec_mass<-read.csv("FAME_molecular_mass_data.csv")
####

FINAL_PLFA<-nmolFAME(data=data,stdRatio=dat_stdRatio,molecular_mass=dat_molec_mass)
write.csv(FINAL_PLFA,file="PLFA_MASTER_nmol_g_soil.csv",row.names=FALSE)

#########################################################################################################
# FUNCTION: getBugGroup
# Appends corresponding microbial group for each FAME compound in a column at the end of the data frame
# input: data = master data frame
# output: Original data frame with microbial group column added to the end
#--------------------------------------------------------------------------------------------------------

getBugGroup<-function(data) {
  N<-length(data[,1])
  group_vec<-rep(NA,times=N)
  
  for (i in 1:N) {
    if(data[i,4]=="12:0"){
      group_vec[i]<-"standard"
    }
    if(data[i,4]=="14:0"){
      group_vec[i]<-"bacteria"
    }
    if(data[i,4]=="a15:0"){
      group_vec[i]<-"gram+bacteria"
    }
    if(data[i,4]=="15:0"){
      group_vec[i]<-"bacteria"
    }
    if(data[i,4]=="16:1n9/i16:0"){
      group_vec[i]<-"bacteria"
    }
    if(data[i,4]=="16:0"){
      group_vec[i]<-"bacteria"
    }
    if(data[i,4]=="i17:0"){
      group_vec[i]<-"gram+bacteria"
    }
    if(data[i,4]=="17:0cy"){
      group_vec[i]<-"gram-bacteria"
    }
    if(data[i,4]=="17:0"){
      group_vec[i]<-"bacteria"
    }
    if(data[i,4]=="18:0"){
      group_vec[i]<-"bacteria"
    }
    if(data[i,4]=="19:0cy"){
      group_vec[i]<-"gram-bacteria"
    }
    if(data[i,4]=="19:0"){
      group_vec[i]<-"standard"
    }
    if(data[i,4]=="20:4n6"){
      group_vec[i]<-"protozoa"
    }
    if(data[i,4]=="20:5n3"){
      group_vec[i]<-"protozoa"
    }
    if(data[i,4]=="18:2n9,12 and cis18:1n9 and 18:3n3"){
      group_vec[i]<-"fungi"
    }
    if(data[i,4]=="18:2n9,12 and cis-18:1n9 and 183n3"){
      group_vec[i]<-"fungi"
    }
    if(data[i,4]=="trans-18:1n9"){
      group_vec[i]<-"gram-bacteria"
    }
    if(data[i,4]=="18:3n6"){
      group_vec[i]<-"fungi"
    }
    if(data[i,4]=="16:1w5"){
      group_vec[i]<-"AMF"
    }
  }
  outDat<-data.frame(data,Microbial_Group=group_vec)
  return(outDat)  
}

#--------------------------------------------------------------------------------------------------------

PLFA_microbe_groups<-getBugGroup(data=FINAL_PLFA)
write.csv(PLFA_microbe_groups,file="MASTER_PLFA_w_MICROBE_GROUP.csv",row.names=FALSE)

#########################################################################################################
# FUNCTION: addMetaData
# Appends corresponding tree type, rhizosphere manipulation treatment, and strach vs leaf substrate added to master data frame
# input: data = master data frame
# output: Original data frame with meta data added in 3 separate columns after sample_name column
#--------------------------------------------------------------------------------------------------------

addMetaData<-function(data) {
  
  N<-length(data[,1])
  treeVec<-rep(NA,times=N)
  rhizVec<-rep(NA,times=N)
  substrateVec<-rep(NA,times=N)

  for (i in 1:N) {
    if(data[i,2]=="01.05.1"){
      treeVec[i]<-"G"
      rhizVec[i]<-"1"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="01.13.1"){
      treeVec[i]<-"P"
      rhizVec[i]<-"1"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="01.13.2"){
      treeVec[i]<-"P"
      rhizVec[i]<-"2"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="01.13.3"){
      treeVec[i]<-"P"
      rhizVec[i]<-"3"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="01.05.3"){
      treeVec[i]<-"G"
      rhizVec[i]<-"3"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="10.02.2"){
      treeVec[i]<-"G"
      rhizVec[i]<-"2"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="10.02.1"){
      treeVec[i]<-"G"
      rhizVec[i]<-"1"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="10.02.3"){
      treeVec[i]<-"G"
      rhizVec[i]<-"3"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="11.06.2"){
      treeVec[i]<-"G"
      rhizVec[i]<-"2"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="11.10.1"){
      treeVec[i]<-"G"
      rhizVec[i]<-"1"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="11.10.3"){
      treeVec[i]<-"G"
      rhizVec[i]<-"3"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="11.06.3"){
      treeVec[i]<-"G"
      rhizVec[i]<-"3"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="12.09.1"){
      treeVec[i]<-"G"
      rhizVec[i]<-"1"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="12.09.2"){
      treeVec[i]<-"G"
      rhizVec[i]<-"2"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="12.01.2"){
      treeVec[i]<-"G"
      rhizVec[i]<-"2"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="12.01.3"){
      treeVec[i]<-"G"
      rhizVec[i]<-"3"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="12.05.2"){
      treeVec[i]<-"G"
      rhizVec[i]<-"2"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="12.05.3"){
      treeVec[i]<-"G"
      rhizVec[i]<-"3"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="12.09.3"){
      treeVec[i]<-"G"
      rhizVec[i]<-"3"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="14.10.3"){
      treeVec[i]<-"G"
      rhizVec[i]<-"3"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="14.10.1"){
      treeVec[i]<-"G"
      rhizVec[i]<-"1"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="14.10.2"){
      treeVec[i]<-"G"
      rhizVec[i]<-"2"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="15.16.2"){
      treeVec[i]<-"P"
      rhizVec[i]<-"2"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="17.01.1"){
      treeVec[i]<-"G"
      rhizVec[i]<-"1"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="17.01.3"){
      treeVec[i]<-"G"
      rhizVec[i]<-"3"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="17.06.1"){
      treeVec[i]<-"P"
      rhizVec[i]<-"1"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="17.06.2"){
      treeVec[i]<-"P"
      rhizVec[i]<-"2"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="17.01.1b"){
      treeVec[i]<-"G"
      rhizVec[i]<-"1"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="02.08.1"){
      treeVec[i]<-"G"
      rhizVec[i]<-"1"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="02.08.2"){
      treeVec[i]<-"G"
      rhizVec[i]<-"2"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="02.08.3"){
      treeVec[i]<-"G"
      rhizVec[i]<-"3"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="02.10.1"){
      treeVec[i]<-"G"
      rhizVec[i]<-"1"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="02.10.3"){
      treeVec[i]<-"G"
      rhizVec[i]<-"3"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="03.07.2"){
      treeVec[i]<-"P"
      rhizVec[i]<-"2"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="03.07.3"){
      treeVec[i]<-"P"
      rhizVec[i]<-"3"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="05.03.1"){
      treeVec[i]<-"P"
      rhizVec[i]<-"1"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="05.03.2"){
      treeVec[i]<-"P"
      rhizVec[i]<-"2"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="05.03.3"){
      treeVec[i]<-"P"
      rhizVec[i]<-"3"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="05.09.2"){
      treeVec[i]<-"P"
      rhizVec[i]<-"2"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="05.09.3"){
      treeVec[i]<-"P"
      rhizVec[i]<-"3"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="06.06.1"){
      treeVec[i]<-"P"
      rhizVec[i]<-"1"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="06.06.2"){
      treeVec[i]<-"P"
      rhizVec[i]<-"2"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="06.06.3"){
      treeVec[i]<-"P"
      rhizVec[i]<-"3"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="17.06.3"){
      treeVec[i]<-"P"
      rhizVec[i]<-"3"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="18.05.1"){
      treeVec[i]<-"P"
      rhizVec[i]<-"1"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="18.05.2"){
      treeVec[i]<-"P"
      rhizVec[i]<-"2"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="18.05.3"){
      treeVec[i]<-"P"
      rhizVec[i]<-"3"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="20.04.1"){
      treeVec[i]<-"P"
      rhizVec[i]<-"1"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="20.04.2"){
      treeVec[i]<-"P"
      rhizVec[i]<-"2"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="22.07.2"){
      treeVec[i]<-"G"
      rhizVec[i]<-"2"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="22.07.3"){
      treeVec[i]<-"G"
      rhizVec[i]<-"3"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="24.02.1"){
      treeVec[i]<-"G"
      rhizVec[i]<-"1"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="24.02.2"){
      treeVec[i]<-"G"
      rhizVec[i]<-"2"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="24.02.3"){
      treeVec[i]<-"G"
      rhizVec[i]<-"3"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="28.01.2"){
      treeVec[i]<-"P"
      rhizVec[i]<-"2"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="28.01.3"){
      treeVec[i]<-"P"
      rhizVec[i]<-"3"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="34.06.1"){
      treeVec[i]<-"G"
      rhizVec[i]<-"1"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="34.06.3"){
      treeVec[i]<-"G"
      rhizVec[i]<-"3"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="37.04.1"){
      treeVec[i]<-"P"
      rhizVec[i]<-"1"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="37.04.2"){
      treeVec[i]<-"P"
      rhizVec[i]<-"2"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="37.04.3"){
      treeVec[i]<-"P"
      rhizVec[i]<-"3"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="43.02.1"){
      treeVec[i]<-"G"
      rhizVec[i]<-"1"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="43.02.2"){
      treeVec[i]<-"G"
      rhizVec[i]<-"2"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="43.02.3"){
      treeVec[i]<-"G"
      rhizVec[i]<-"3"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="46.03.2"){
      treeVec[i]<-"P"
      rhizVec[i]<-"2"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="46.03.3"){
      treeVec[i]<-"P"
      rhizVec[i]<-"3"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="48.01.1"){
      treeVec[i]<-"G"
      rhizVec[i]<-"1"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="48.01.3"){
      treeVec[i]<-"G"
      rhizVec[i]<-"3"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="52.04.1"){
      treeVec[i]<-"P"
      rhizVec[i]<-"1"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="52.04.2"){
      treeVec[i]<-"P"
      rhizVec[i]<-"2"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="52.04.3"){
      treeVec[i]<-"P"
      rhizVec[i]<-"3"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="55.11.1"){
      treeVec[i]<-"P"
      rhizVec[i]<-"1"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="55.11.2"){
      treeVec[i]<-"P"
      rhizVec[i]<-"2"
      substrateVec[i]<-"S"
    }
    if(data[i,2]=="58.05.1"){
      treeVec[i]<-"P"
      rhizVec[i]<-"1"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="58.05.2"){
      treeVec[i]<-"P"
      rhizVec[i]<-"2"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="58.05.3"){
      treeVec[i]<-"P"
      rhizVec[i]<-"3"
      substrateVec[i]<-"L"
    }
    if(data[i,2]=="63.04.2"){
      treeVec[i]<-"P"
      rhizVec[i]<-"2"
      substrateVec[i]<-"L"
    } # close final if statement
  } # for loop
  rhizVec<-as.factor(rhizVec)
  outDat<-data.frame(ID=data[,1],SampleName=data[,2],Tree_Type=treeVec,Rhizosphere_Manipulation=rhizVec,Substrate_Type=substrateVec,data[,3:12])
  return(outDat)
}

#--------------------------------------------------------------------------------------------------------

PLFA_bugs<-read.csv("MASTER_PLFA_w_MICROBE_GROUP.csv")
PLFA_with_meta<-addMetaData(data=PLFA_bugs)
write.csv(PLFA_with_meta,file="FINAL_PLFA_with_metadata.csv",row.names=FALSE)

