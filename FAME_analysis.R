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

#########################################################################################################
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

# percentC()

# fame_molec<-read.csv("fame_percent_C.csv")
# fameCog<-percentC(data=fame_molec) # Data frame with percent_C calculations in last column
# fameC<-data.frame(fameCog[,1:5],percent_C=fameCog[,7]) # Remove weird blank column
# write.csv(fameC,file="FAME_carbon.csv")

FAMEcarbon<-read.csv("FAME_carbon.csv")
carbPLFA<-data.frame(FAME_compound=FAMEcarbon[,2],FAME_per_C=FAMEcarbon[,7])

##### Creating dataframe that includes only variable relevant for analysis

cleanDat<-read.csv("CleanFAME_dat.csv") # This dataset has been edited using RegEx in BBedit to make compound labels consistent (ex. 14:00 and 14:0 = 14:0) AND sample names consitent (format = ##.##.#)

datPLFA<-data.frame(ID=seq(from=1,to=length(cleanDat$X)),cleanDat[,2:4],cleanDat[7],cleanDat$Area.All,cleanDat$d.13C.12C,cleanDat$AT..13C.12C)

names(datPLFA)<-c("ID","SampleName","InjectionVol_ul","FAME_compound","Sample_Dilution","Peak_Area","d13","Atom_per_13C")

#########################################################################################################
# FUNCTION: includeCarbon
# Append FAME percent C to the PLFA data frame
# input: perC = 2-column dataframe with FAME compound name and percent carbon 0.00)
# dat = dataframe that percent C will be appended to where first column= ID
# output: Data frame with percent C of each compound in last column
#--------------------------------------------------------------------------------------------------------

includeCarbon<-function(perC,data) {
  N<-length(data[,1])
  perCvect<-rep("NA",length=N)
  
  for(i in 1:N) {
    if(data[i,4]=="12:0"){
      perCvect[i]<-round(perC[1,2],digits=4)
    }
    if(data[i,4]=="14:0"){
      perCvect[i]<-round(perC[2,2],digits=4)
    }
    if(data[i,4]=="a15:0"){
      perCvect[i]<-round(perC[3,2],digits=4)
    }
    if(data[i,4]=="15:0"){
      perCvect[i]<-round(perC[4,2],digits=4)
    }
    if(data[i,4]=="16:1n9/i16:0"){
      perCvect[i]<-round(perC[5,2],digits=4)
    }
    if(data[i,4]=="16:0"){
      perCvect[i]<-round(perC[6,2],digits=4)
    }
    if(data[i,4]=="i17:0"){
      perCvect[i]<-round(perC[7,2],digits=4)
    }
  
  }
  outDat<-data.frame(data,FAME_per_C=perCvect)
  return(head(outDat))
}

testerDAT<-datPLFA
testerCARB<-carbPLFA

includeCarbon(perC=testerCARB,data=testerDAT)
