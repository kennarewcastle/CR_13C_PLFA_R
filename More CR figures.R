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


# Root biomass by rhizosphere treatment ----------------------------------
roots_rhizo<-lm(data$Dry_root_wt~data$Exclusion)
Anova(roots_rhizo)
ggplot(data,aes(x=Exclusion,y=Dry_root_wt,fill=Tree_species)) +
  geom_boxplot(lwd=1.5) +
  ylab(expression(bold(paste("Root Biomass (dry mg)")))) +
  xlab(label=expression(bold("Rhizosphere Manipulation"))) +
  labs(fill="Tree Species") +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme_classic()

ggplot(data=data, aes(data$Dry_root_wt,fill=Exclusion)) + 
  geom_density(alpha=0.7) +
  labs(fill="Exclusion",x="Dry Root Weight (mg)",y="Density") +
  scale_fill_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme_classic()

# C-enzymes by rhizosphere manipulation -----------------------------------

# ANOVA for C-enzymes by tree species
tree_Cenzy<-lm(data$C_enzymes~data$Tree_species)
Anova(tree_Cenzy) # No difference in C-enzyme activity between tree species
ggplot(data,aes(x=Tree_species,y=C_enzymes)) +
  geom_boxplot(lwd=1.3) +
  ylab(label="Carbon Acquiring Enzyme Activity") +
  xlab(label="Tree Species") +
  theme_classic()

# ANOVA for C-enzymes by rhizosphere manipulation
rhizo_Cenzy<-lm(data$C_enzymes~data$Exclusion)
Anova(rhizo_Cenzy) # No difference in C-enzyme activity between rhizosphere manipulation

Cenzy<-ggplot(data,aes(x=Exclusion,y=C_enzymes,fill=Tree_species)) +
  geom_boxplot(lwd=1.3) +
  ylab(expression(bold(paste("C-Acquiring Enzyme Activity"," (nmol"," g"^-1," h"^-1,")")))) +
  xlab(label=expression(bold("Rhizosphere Manipulation"))) +
  labs(fill="Tree Species") +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme_classic()

ggsave(filename="Cenzy_rhizo.jpg")

rhizo_Cenzy_reg<-lm(data$C_enzymes~data$Dry_root_wt)
summary(rhizo_Cenzy_reg)
ggplot(data,aes(x=Dry_root_wt,y=C_enzymes)) +
  geom_point()+
  stat_smooth(method="lm",formula=y~x,fullrange=TRUE,colour="black",fill="black",size=1.25) +
  ylab(expression(bold(paste("C-Acquiring Enzyme Activity"," (nmol"," g"^-1," h"^-1,")")))) +
  xlab(label=expression(bold("Root Dry Weight (mg)"))) +
  theme_classic()


# Nutrient enzymes by rhizosphere manipulation ----------------------------

# ANOVA for nut-enzymes by tree species
tree_nutenzy<-lm(data$nut_enzymes~data$Tree_species)
Anova(tree_nutenzy) # No difference in nut-enzyme activity between tree species
ggplot(data,aes(x=Tree_species,y=nut_enzymes)) +
  geom_boxplot(lwd=1.3) +
  ylab(label="Nutrient Acquiring Enzyme Activity") +
  xlab(label="Tree Species") +
  theme_classic()

# ANOVA for C-enzymes by rhizosphere manipulation
rhizo_nutenzy<-lm(data$nut_enzymes~data$Exclusion)
Anova(rhizo_nutenzy) # No difference in C-enzyme activity between rhizosphere manipulation
nutenzy<-ggplot(data,aes(x=Exclusion,y=nut_enzymes,fill=Tree_species)) +
  geom_boxplot(lwd=1.3) +
  ylab(expression(bold(paste("Nutrient-Acquiring Enzyme Activity"," (nmol"," g"^-1," h"^-1,")")))) +
  xlab(label=expression(bold("Rhizosphere Manipulation"))) +
  guides(fill=FALSE) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme_classic()

#ggsave(filename="nutenzy_rhizo.jpg")

grid.arrange(nutenzy,Cenzy,nrow=1)
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

# Data separation for starch and leaf labeled cores
starch_dat<-filter(data,Isotope_label=="S")
leaf_dat<-filter(data,Isotope_label=="L")

# ANOVA for d13 by tree species for starch
tree_starch_co2<-lm(starch_dat$d13_d5~starch_dat$Tree_species)
Anova(tree_starch_co2) # No difference in starch CO2 d13 between tree species (p=0.7339)
d13_tree_starch<-ggplot(starch_dat,aes(x=Tree_species,y=d13_d5)) +
  geom_boxplot() +
  ylab(label=expression(bold(paste("\u03B4"^{bold("13")}, "CO"[bold("2")]," (\u2030)")))) +
  xlab(label=expression(bold("Tree Species"))) +
  scale_x_discrete(labels=c("Goeth"="Goethalsia","Penta"="Pentaclethra")) +
  ylim(-25,50) +
  labs(title = expression(bold("Starch Substrate"))) +
  theme_classic() +
  theme(axis.text=element_text(colour="black",size=12))

# ANOVA for d13 by tree species for leaf
tree_leaf_co2<-lm(leaf_dat$d13_d9~leaf_dat$Tree_species)
Anova(tree_leaf_co2) # Significant difference in d13 in CO2 by tree species for leaf (p=0.03057)
d13_tree_leaf<-ggplot(leaf_dat,aes(x=Tree_species,y=d13_d9)) +
  geom_boxplot() +
  ylab(label=expression(bold(paste("\u03B4"^{bold("13")}, "CO"[bold("2")]," (\u2030)")))) +
  xlab(label=expression(bold("Tree Species"))) +
  scale_x_discrete(labels=c("Goeth"="Goethalsia","Penta"="Pentaclethra")) +
  labs(title = expression(bold("Leaf Substrate"))) +
  annotate("text", x = 2, y = 92, label = "*p = 0.031", size=4) +
  theme_classic() +
  theme(axis.text=element_text(colour="black",size=12))

# Two paneled figure for effect of tree species on 13C-CO2 for starch and leaf substrates
d13_tree<-grid.arrange(d13_tree_starch,d13_tree_leaf,nrow=1)
ggsave(filename="d13_CO2_tree.jpg",plot=d13_tree)



###########################
# ANOVA for d13 by rhizosphere manipulation for 
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


# Data read in--> starch and leaf will be separated
d13PLFAdat<-read.csv("FINAL_PLFA_with_metadata.csv")
d13PLFAdat$Rhizosphere_Manipulation<-as.factor(d13PLFAdat$Rhizosphere_Manipulation)
PLFA_starch<-filter(d13PLFAdat,Substrate_Type=="S")
PLFA_leaf<-filter(d13PLFAdat,Substrate_Type=="L")

