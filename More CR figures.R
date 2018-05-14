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


# C-enzymes by rhizosphere manipulation -----------------------------------

data<-read.csv("Costa Rica Master sheet annotated.csv")
data$Exclusion<-as.factor(data$Exclusion) # Rhizosphere manipulation is now a factor instead of numeric
data<-filter(data,Exclusion!="NA")
C_enzymes<-data$cbh + data$a_gluc + data$b_gluc
nut_enzymes<-data$nag + data$phos + data$lap
data<-data.frame(data,C_enzymes,nut_enzymes)

# ANOVA for tree_species by C-enzymes
tree_Cenzy<-lm(data$C_enzymes~data$Tree_species)
Anova(tree_Cenzy) # No difference in C-enzyme activity between tree species
ggplot(data,aes(x=Tree_species,y=C_enzymes)) +
  geom_boxplot() +
  ylab(label="Carbon Degrading Enzyme Activity") +
  xlab(label="Tree Species") +
  theme_classic()
