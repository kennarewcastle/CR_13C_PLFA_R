### Calculation of carbon use efficiency from 13CO2 and 13C PLFA
### KER
### Created on: 16 November 2021
### Last modified:
###         30 November 2021: Summed 13C-PLFA mass across all FAME compounds per sample
###         30 November 2021: Added analysis for total 13C-CO2 for each core
# Load packages -----------------------------------------------------------

library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

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
  
  outDat1<-data.frame(data,"StandardRatio"=stdVec)
  
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

data$InjectionVol_ul[1239:1255] <- "3 ul" # Fixes weird typo in convention
data<-data %>% separate(InjectionVol_ul,into=c("Injection_vol_ul",NA),sep=" ",remove=TRUE,convert=FALSE,extra="merge",fill="warn") # Pulls out the number from the 3 ul character string
data$Injection_vol_ul<-as.numeric(data$Injection_vol_ul)

stdRatio<-read.csv("PeakAreaC_ratio_per_sample.csv")

C_data<-ug_C_FAME(data=data,stdRatio=dat_stdRatio)

# write.csv(C_data,file="Total_g_FAME_for_specific_resp.csv",row.names=FALSE)

# Multiply atom % by C content to get mass labeled C ----------------------

C_data<-read.csv("Total_g_FAME_for_specific_resp.csv")
C_data<-data.frame(C_data$ID,C_data$SampleName,C_data$FAME_compound,C_data$d13,C_data$Atom_per_13C,C_data$Total_C_ug)
names(C_data)<-c("Unique_ID","Core_ID","FAME_Compound","d13","Atom_per_13C","Total_C_ug")

C_data$Mass_13C_ug<-(C_data$Atom_per_13C/100)*C_data$Total_C_ug

# Aggregate 13C per core across all FAME compounds ------------------------

outdat<-ddply(C_data, .(Core_ID), summarise, Total_13C_ug=sum(Mass_13C_ug))
write.csv(outdat,file="Total_13C_per_Core.csv",row.names=FALSE) # Core leaf/starch label type added manually after writing this file


# Figures for microbial 13C by core ---------------------------------------

data<-read.csv("Total_13C_per_Core.csv")
data$Mesh_Type<-as.factor(data$Mesh_Type)

leaf<-filter(data,Label_Type=="L")
starch<-filter(data,Label_Type=="S")

### These ylim scalers eliminate 1 outlier in starch with total ug 13C > 3000, change for publication
leaf_rhizo_13C<-ggplot(leaf,aes(x=Mesh_Type,y=Total_13C_ug)) +
  geom_boxplot(lwd=1.5,outlier.shape=NA) +
  geom_point(size=5,alpha=0.7) +
  ylim(190,2500) +
  labs(y=expression(bold(paste("Total Microbial PLFA"," "^13,"C (",mu,"g)"))),x="Rhizosphere Exclusion",title="Leaf") +
  scale_x_discrete(labels=c("-R-M","-R+M","+R+M")) +
  theme_classic()+
  theme(
    axis.text.x=element_text(colour="black",size=14),
    axis.text.y=element_text(colour="black"),
    axis.title.x=element_text(face="bold",size=16),
    axis.title.y=element_text(face="bold",size=16),
    plot.title=element_text(face="bold",size=18))

starch_rhizo_13C<-ggplot(starch,aes(x=Mesh_Type,y=Total_13C_ug)) +
  geom_boxplot(lwd=1.5,outlier.shape=NA) +
  geom_point(size=5,alpha=0.7) +
  ylim(190,2500) +
  labs(y=expression(bold(paste("Total Microbial PLFA"," "^13,"C (",mu,"g)"))),x="Rhizosphere Exclusion",title="Starch") +
  scale_x_discrete(labels=c("-R-M","-R+M","+R+M")) +
  theme_classic()+
  theme(
    axis.text.x=element_text(colour="black",size=14),
    axis.text.y=element_text(colour="black"),
    axis.title.x=element_text(face="bold",size=16),
    axis.title.y=element_text(face="bold",size=16),
    plot.title=element_text(face="bold",size=18))

grid.arrange(starch_rhizo_13C,leaf_rhizo_13C,ncol=2)


# Calculate mass 13C-CO2 per core -----------------------------------------

### This conversion is kind of tricky... First use the ideal gas law with the following variables to solve for the number of moles of gas total in each core (volume changes for each core with the depth of the headspace in teh cores). P = atmospheric pressure at 35 m above sea level (from La Selva website) = 345.39 Pa, T= 298.91 K (average of all temperatures taken inside each core in the field notes spreadsheet), V = units are m3, calculated by multiplying the area of the cores given the cores' radii (2.5 cm) by the depth of soil from the top of the cores (field notes spreadsheet), R = gas constant (8.314510).

CO2<-read.csv("CR2016_13CO2_Gas_Samples.csv")

# Calculate gas volume
CO2$Core_Volume_m3<-((pi*(2.5)^2)*CO2$Core_Depth_cm)/1000000

# Calculate moles of gas = PV/RT
CO2$Gas_Moles<-(345.39*CO2$Core_Volume_m3)/(8.314510*298.91)

# Multiply moles of gas by ppmv/1000000 to get moles of CO2 on each day, then multiply by 13 (molecular weight of 13C carbon) to get g of 13C emitted by each core over 10 minute capping interval, then multiply by 1000000 to get ug 13C
CO2$C13_ug_D1<-CO2$Gas_Moles*(CO2$ppmv_day1/1000000)*100000
CO2$C13_ug_D2<-CO2$Gas_Moles*(CO2$ppmv_day2/1000000)*100000
CO2$C13_ug_D3<-CO2$Gas_Moles*(CO2$ppmv_day3/1000000)*100000
CO2$C13_ug_D4<-CO2$Gas_Moles*(CO2$ppmv_day4/1000000)*100000
CO2$C13_ug_D5<-CO2$Gas_Moles*(CO2$ppmv_day5/1000000)*100000
CO2$C13_ug_D6<-CO2$Gas_Moles*(CO2$ppmv_day6/1000000)*100000
CO2$C13_ug_D7<-CO2$Gas_Moles*(CO2$ppmv_day7/1000000)*100000
CO2$C13_ug_D9<-CO2$Gas_Moles*(CO2$ppmv_day9/1000000)*100000

# Separate leaf data from starch data because cores were sampled/harvested on different dates
leaf<-filter(CO2,Label_Type=="leaf") # Gas sampling days are 1, 4, 6, 7, 9
starch<-filter(CO2,Label_Type=="starch") # Gas sampling days are 1, 2, 3, 4, 5

# Check to see if any starch cores are missing gas samples for certain days... Missing samples will be calculated manually after calculating total mass produced for the rest of the dataset

length(starch$ppmv_day1[is.na(starch$ppmv_day1)]) # All samples present day 1

length(starch$ppmv_day2[is.na(starch$ppmv_day2)]) # Missing 2 samples (6 and 15)
  starch[c(6,15),] # Missing samples are 10.2.3 and 17.1.3

length(starch$ppmv_day3[is.na(starch$ppmv_day3)]) # All samples present day 3

length(starch$ppmv_day4[is.na(starch$ppmv_day4)]) # Missing 1 sample (3)
  starch[3,] # Missing sample is 1.13.3
  
length(starch$ppmv_day5[is.na(starch$ppmv_day5)]) # Missing 1 sample (45)
  starch[45,] # Missing sample is 52.4.3
  
# Calculate total 13C emitted over 5 day starch incubation (considering 10 minute cap period at time of gas sample collection)... 1440 = number of minutes in 24 hours
starch$Total_13C_ug<-(starch$C13_ug_D1/10*1440) + (starch$C13_ug_D2/10*1440) + (starch$C13_ug_D3/10*1440) + (starch$C13_ug_D4/10*1440) + (starch$C13_ug_D5/10*1440)
  
# Redo cores with missing samples manually
  
  # Cores with day 2 samples missing (day 2 ppmv = average ug of day 1 and day 3)
  starch$Total_13C_ug[6]<-(starch$C13_ug_D1[6]/10*1440) + (mean(c(starch$C13_ug_D1[6],starch$C13_ug_D3[6]))/10*1440) + (starch$C13_ug_D3[6]/10*1440) + (starch$C13_ug_D4[6]/10*1440) + (starch$C13_ug_D5[6]/10*1440)
  
  starch$Total_13C_ug[15]<-(starch$C13_ug_D1[15]/10*1440) + (mean(c(starch$C13_ug_D1[15],starch$C13_ug_D3[15]))/10*1440) + (starch$C13_ug_D3[15]/10*1440) + (starch$C13_ug_D4[15]/10*1440) + (starch$C13_ug_D5[15]/10*1440)
  
  # Core with day 4 sample missing
  starch$Total_13C_ug[3]<-(starch$C13_ug_D1[3]/10*1440) + (starch$C13_ug_D2[3]/10*1440) + (starch$C13_ug_D3[3]/10*1440) + (mean(c(starch$C13_ug_D1[15],starch$C13_ug_D3[15]))/10*1440) + (starch$C13_ug_D5[3]/10*1440)

  


