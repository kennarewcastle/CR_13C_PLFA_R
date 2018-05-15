# More visualizations for Costa Rica Data
# May 14 2018
# KER


# Preliminaries -----------------------------------------------------------

library(dplyr)
library(ggplot2)
library(gridExtra)
library(car)

# Trying to use endophyte colonization data as a proxy for mycorrhizal biomass ------------------------

endoDat<-read.csv("CR_2016_fungal_colonization_DOE_cores.csv")
names(endoDat)
endoDat$Core<-as.factor(endoDat$Core)
endoDat<-filter(endoDat,Core!="NA")

ggplot(endoDat,aes(x=Core,y=per_AMF)) +
  geom_boxplot()

mod<-lm(endoDat$per_AMF~endoDat$Core)
anova(mod)
anovaMod<-aov(mod)
TukeyHSD(anovaMod)

names(endoDat)

DatDat<-read.csv("Costa Rica Master sheet annotated.csv")
names(DatDat)
hyphal_scaled<-DatDat$Dry_root_wt*DatDat$per_AMF
newDat<-data.frame(DatDat,hyphal_scaled)
newDat$Exclusion<-as.factor(newDat$Exclusion)
newDat<-filter(newDat,Exclusion!="NA")

ggplot(newDat,aes(x=Exclusion,y=hyphal_scaled)) + # y-axis is the mass of endophyte-colonized roots in each core
  geom_boxplot()


# Data for summary figures (ALL CORES) ------------------------------------------------
data<-read.csv("Costa Rica Master sheet annotated.csv")
data$Exclusion<-as.factor(data$Exclusion) # Rhizosphere manipulation is now a factor instead of numeric
data<-filter(data,Exclusion!="NA")
C_enzymes<-data$cbh + data$a_gluc + data$b_gluc # Creating binned C and nutrient enzyme variables
nut_enzymes<-data$nag + data$phos + data$lap
data<-data.frame(data,C_enzymes,nut_enzymes)

# C-enzymes by rhizosphere manipulation -----------------------------------

# ANOVA for C-enzymes by tree species
tree_Cenzy<-lm(data$C_enzymes~data$Tree_species)
Anova(tree_Cenzy) # No difference in C-enzyme activity between tree species
ggplot(data,aes(x=Tree_species,y=C_enzymes)) +
  geom_boxplot() +
  ylab(label="Carbon Acquiring Enzyme Activity") +
  xlab(label="Tree Species") +
  theme_classic()

# ANOVA for C-enzymes by rhizosphere manipulation
rhizo_Cenzy<-lm(data$C_enzymes~data$Exclusion)
Anova(rhizo_Cenzy) # No difference in C-enzyme activity between rhizosphere manipulation
ggplot(data,aes(x=Exclusion,y=C_enzymes)) +
  geom_boxplot(lwd=1.3) +
  ylab(expression(bold(paste("C-Acquiring Enzyme Activity"," (nmol"," g"^-1," h"^-1,")")))) +
  xlab(label=expression(bold("Rhizosphere Manipulation"))) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme_classic()

ggsave(filename="Cenzy_rhizo.jpg")


# Nutrient enzymes by rhizosphere manipulation ----------------------------

# ANOVA for nut-enzymes by tree species
tree_nutenzy<-lm(data$nut_enzymes~data$Tree_species)
Anova(tree_nutenzy) # No difference in nut-enzyme activity between tree species
ggplot(data,aes(x=Tree_species,y=nut_enzymes)) +
  geom_boxplot() +
  ylab(label="Nutrient Acquiring Enzyme Activity") +
  xlab(label="Tree Species") +
  theme_classic()

# ANOVA for C-enzymes by rhizosphere manipulation
rhizo_nutenzy<-lm(data$nut_enzymes~data$Exclusion)
Anova(rhizo_nutenzy) # No difference in C-enzyme activity between rhizosphere manipulation
ggplot(data,aes(x=Exclusion,y=nut_enzymes)) +
  geom_boxplot(lwd=1.3) +
  ylab(expression(bold(paste("Nutrient-Acquiring Enzyme Activity"," (nmol"," g"^-1," h"^-1,")")))) +
  xlab(label=expression(bold("Rhizosphere Manipulation"))) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme_classic()

ggsave(filename="nutenzy_rhizo.jpg")


# Respiration by rhizosphere manipulation ---------------------------------

# Respiration data on harvest day is different for starch and leaf cores, create new column for harvest day resp that contains appropriate resp measurement for each core (starch = d5, leaf = d9)

#########################################################################################################
# FUNCTION: HarvestResp
# Creates vector of respiration value on the day each core was harvested, appends to dataframe
# input: data= dataframe
# output: dataframe with harvestresp column appended to end
#--------------------------------------------------------------------------------------------------------

HarvestResp<-function(data) {
  N<-length(data$Exclusion)
  harVec<-rep(NA,times=N)

  for (i in 1:N) {
    if (data[i,4]=="L") {
      harVec[i]<-data$Efflux_d9[i]
    }
    if (data[i,4]=="LW") {
      harVec[i]<-data$Efflux_d9[i]
    }
    if (data[i,4]=="S") {
      harVec[i]<-data$Efflux_d5[i]
    }
    if (data[i,4]=="SW") {
      harVec[i]<-data$Efflux_d5[i]
    }
  }
  outDat<-data.frame(data,harvest_day_resp=harVec)
  return(outDat)
}
#--------------------------------------------------------------------------------------------------------

data<-HarvestResp(data=data) # Respiration rate on harvest day is now in "data" dataframe

# ANOVA for respiration by tree species
tree_resp<-lm(data$harvest_day_resp~data$Tree_species)
Anova(tree_resp) # No difference in resp between tree species (p=0.2191)
ggplot(data,aes(x=Tree_species,y=nut_enzymes)) +
  geom_boxplot() +
  ylab(label="Respiration Rate") +
  xlab(label="Tree Species") +
  theme_classic()

# ANOVA for C-enzymes by rhizosphere manipulation
rhizo_resp<-lm(data$harvest_day_resp~data$Exclusion)
Anova(rhizo_resp) # Marginal differences in resp between rhizosphere treatments (p=0.05852)
aov_rhizo_resp<-aov(rhizo_resp)
TukeyHSD(aov_rhizo_resp) # 1 is must diferent from 3
ggplot(data,aes(x=Exclusion,y=harvest_day_resp)) +
geom_boxplot(lwd=1.3) +
  ylab(expression(bold(paste("Respiration (",mu,"mol"," ","CO"[2]," ","m"^-2," s"^-1,")")))) +
  xlab(label=expression(bold("Rhizosphere Manipulation"))) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme_classic()

ggsave(filename="resp_rhizo.jpg")


# d13 by rhizosphere manipulation for leaf and starch ---------------------




# Data read in--> starch and leaf will be separated
d13dat<-read.csv("FINAL_PLFA_with_metadata.csv")
d13dat$Rhizosphere_Manipulation<-as.factor(d13dat$Rhizosphere_Manipulation)
starch_dat<-filter(d13dat,Substrate_Type=="S")
leaf_dat<-filter(d13dat,Substrate_Type=="L")

# ANOVA 
