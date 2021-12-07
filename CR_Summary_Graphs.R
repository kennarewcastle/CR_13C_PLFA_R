# Costa Rica Summary Graphs
# Created: 12 December 2018
# Modified: 
#     14 January 2019: Summary graphs added.
#     08 April 2021: Revisiting this script, making sure code runs appropriately.
# KER


# Packages ----------------------------------------------------------------
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(reshape2)
library(stringr)

# Global variables --------------------------------------------------------
data<-read.csv("FINAL_PLFA_with_metadata.csv")
data<-data[,2:15] # Gets rid of ID column
data<-filter(data,Microbial_Group!="standard")
data$Rhizosphere_Manipulation<-as.factor(data$Rhizosphere_Manipulation)
data$umol_FAME<-data$nmol_FAME_per_g_soil/1000

master<-read.csv("Costa Rica Master sheet annotated.csv")


# Append PLFA data to master ----------------------------------------------

# Sort data frame so compounds are listed in same order, NA used for samples for which certain compounds weren't detected.
data<-data[order(data$SampleName,data$FAME_compound),]

# Remove standards from data frame
data<-filter(data,"Microbial_Group"!="standard")

# FAME CONCENTRATION --- Strip "data" dataframe down to sample name, compound name, compound concentration to reformat in wide form
data1<-data.frame("Sample_Name"=data[,1],"FAME_Compound"=data[,6],"umol_FAME"=data[,15])

# Use dcast in the reshape2 package to convert long to wide for dataframes

# This function checks the number of concentration values for each compound to check for duplicates.
data1_wide<-dcast(data1,Sample_Name ~ FAME_Compound, value.var="umol_FAME", fun.aggregate = length)

# Looks like there are duplicates for most of the compounds for 11.06.2 and 11.06.3, but they're identical duplicates, so using the mean to aggregate the compound concentrations for these samples will eliminate the duplication issue. Using the 'fill = 0' argument puts zero for concentration of undetected compounds instead of Na in order to calculate mean.
data1_wide<-dcast(data1,Sample_Name ~ FAME_Compound, value.var="umol_FAME", fun.aggregate = mean, fill=0)

# Calculate saturated/(unsaturated + cyclopropyl) 
# In fatty acid terminology, using Guckert et al 29184 appl. environ. microbiol.
.<-data1_wide
data1_wide$Saturated<- .[,2] + .[,3] + .[,4] + .[,7] + .[,9]
data1_wide$Unsaturated<- .[,5] + .[,6] + .[,8] + .[,10] + .[,11] + .[,12] + .[,13] + .[,14] + .[,15] + .[,16] + .[,17]
data1_wide$Stress_Ratio<-data1_wide$Saturated/data1_wide$Unsaturated
data1_wide$Rhizo_Manipulation<-str_sub(data1_wide$Sample_Name,-1,-1)

# FAME d13 --- Strip 'data' down to samplename, compound name, d13
data2<-data.frame("Sample_Name"=data[,1],"FAME_Compound"=data[,6],"FAME_d13"=data[,9])

# Reshape from long to wide format.
data2_wide<-dcast(data2,Sample_Name ~ FAME_Compound, value.var="FAME_d13", fun.aggregate=mean, na.rm=TRUE)

# Convert NaNs (not a real number from using mean function on NA values) to NA.
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

data2_wide[is.nan(data2_wide)] <- NA

# CONCENTRATION BY MICROBIAL GROUP --- list concentrations of FAME aggregated across microbial group, use this data to compute fungi:bacteria ratios.

data3<-data.frame("Sample_Name"=data[,1],"Microbial_Group"=data[,14],"umol_FAME"=data[,15])

# Reshape from long to wide format, sum concentrations across compounds linked to the same microbial groups
data3_wide<-dcast(data3,Sample_Name ~ Microbial_Group, value.var="umol_FAME", fun.aggregate=sum, fill=0)

.<-data3_wide
fungi_bact<-((.[,2] + .[,4])/(.[,3] + .[,5] + .[,6])) # (AMF + fungi) / (bacteria + g+bact + g-bact)

data3_wide<-data.frame(data3_wide, "fungi_bacteria" = fungi_bact)

colnames(data3_wide)<-c("Sample_Name","AMF","bacteria","fungi","gram_neg_bacteria","gram_pos_bacteria","protozoa","fungi_bacteria_ratio") 

# D13 BY MICROBIAL GROUP --- list d13 of FAME compounds averaged across all compounds in a microbial group, use this data to compute ACTIVE fungi:bacteria ratios.

data4<-data.frame("Sample_Name"=data[,1],"Microbial_Group"=data[,14],"d13"=data[,9])

# Reshape from long to wide format, average d13 values across compounds linked to the same microbial groups
data4_wide<-dcast(data4,Sample_Name ~ Microbial_Group, value.var="d13", fun.aggregate=mean, na.rm=TRUE)

# Convert NaNs (not a real number from using mean function on NA values) to NA.
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

data4_wide[is.nan(data4_wide)] <- NA

.<-data4_wide

d13_fungi_bact<-(((.[,2]+100) + (.[,4]+100))/((.[,3]+100) + (.[,5]+100) + (.[,6]+100))) # (AMF + fungi) / (bacteria + g+bact + g-bact), all values scaled by adding 100 to remove negative values.

data4_wide<-data.frame(data4_wide, "d13_fungi_bacteria" = d13_fungi_bact)

colnames(data4_wide)<-c("Sample_Name","d13_AMF","d13_bacteria","d13_fungi","d13_gram_neg_bacteria","d13_gram_pos_bacteria","d13_protozoa","d13_fungi_bacteria_ratio") 

# Microbial community by rhizosphere and tree -----------------------------
# Two panneled figure, one for each tree type.

data_goeth<-filter(data,Tree_Type=="G")
FAME_rhizo_g<-ggplot(data_goeth,aes(x=Rhizosphere_Manipulation,y=umol_FAME,fill=Microbial_Group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  ylab(label="FAME concentration (\u03BCmol FAME per g dry soil)") +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  xlab(label="Rhizosphere Manipulation") +
  ggtitle("Goethalsia") +
  labs(fill="Microbial Group") +
  ylim(0,42) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))


data_pent<-filter(data,Tree_Type=="P")
FAME_rhizo_p<-ggplot(data_pent,aes(x=Rhizosphere_Manipulation,y=umol_FAME,fill=Microbial_Group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  ylab(label="FAME concentration (\u03BCmol FAME per g dry soil)") +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  xlab(label="Rhizosphere Manipulation") +
  ggtitle("Pentaclethera") +
  labs(fill="Microbial Group") +
  ylim(0,42)+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

grid.arrange(FAME_rhizo_p,FAME_rhizo_g,nrow=1)


# Active microbial community (d13 PLFA) by rhizosphere and tree -----------
# Four panneled figure, one for each substrate x tree species

data_goeth_starch<-filter(data_goeth,Substrate_Type=="S")
d13_g_starch<-ggplot(data_goeth_starch,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_brewer(palette = "Set1") +
  ylab(label=expression(bold(paste("\u03B4"^{bold("13")}, "C PLFA"," (\u2030)")))) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  xlab(label="Rhizosphere Manipulation") +
  ggtitle("Goethalsia | Starch Substrate") +
  labs(fill="Microbial Group") +
  ylim(-75,450) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA),
        legend.position="none")

data_goeth_leaf<-filter(data_goeth,Substrate_Type=="L")
d13_g_leaf<-ggplot(data_goeth_leaf,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_brewer(palette = "Set1") +
  ylab(label=expression(bold(paste("\u03B4"^{bold("13")}, "C PLFA"," (\u2030)")))) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  xlab(label="Rhizosphere Manipulation") +
  ggtitle("Goethalsia | Leaf Substrate") +
  labs(fill="Microbial Group") +
  ylim(-75,450) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA),
        legend.position="none")


data_pent_starch<-filter(data_pent,Substrate_Type=="S")
d13_p_starch<-ggplot(data_pent_starch,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_brewer(palette = "Set1") +
  ylab(label=expression(bold(paste("\u03B4"^{bold("13")}, "C PLFA"," (\u2030)")))) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  xlab(label="Rhizosphere Manipulation") +
  ggtitle("Pentaclethera | Starch Substrate") +
  labs(fill="Microbial Group") +
  ylim(-75,450) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA),
        legend.position="none")

data_pent_leaf<-filter(data_pent,Substrate_Type=="L")
d13_p_leaf<-ggplot(data_pent_leaf,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_brewer(palette = "Set1") +
  ylab(label=expression(bold(paste("\u03B4"^{bold("13")}, "C PLFA"," (\u2030)")))) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  xlab(label="Rhizosphere Manipulation") +
  ggtitle("Pentaclethera | Leaf Substrate") +
  labs(fill="Microbial Group") +
  ylim(-75,450) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA),
        legend.position="none")

grid.arrange(d13_g_starch,d13_g_leaf,d13_p_starch,d13_p_leaf,nrow=2)


# Mycorrhizal abundance by rhizosphere manipulation -----------------------

AM_data<-filter(data,Microbial_Group=="AMF")
AM_rhizo<-ggplot(AM_data,aes(x=Rhizosphere_Manipulation,y=umol_FAME,fill=Tree_Type)) +
  geom_boxplot(lwd=1.5) +
  ylab(label="AMF biomass (\u03BCmol FAME per g dry soil)") +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  xlab(label="Rhizosphere Manipulation") +
  labs(fill="Tree Species") +
  scale_fill_discrete(labels=c("G"="Goethalsia","P"="Pentaclethra")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

myco_mod<-lm(AM_data$umol_FAME~AM_data$Rhizosphere_Manipulation*AM_data$Tree_Type)
anova(myco_mod)

# Dry root mass by DBH -----------------------------------------------------
master3<-filter(master,Exclusion=="3")
DBH_roots_fig<-ggplot(master3,aes(x=Tree_DBH,y=Dry_root_wt)) +
  geom_point() +
  geom_smooth(method=lm) +
  #ylab(label="DBH (cm)") +
  #xlab(label="Tree Species") +
  #scale_x_discrete(labels=c("Goeth"="Goethalsia","Penta"="Pentaclethra"),breaks=c("Goeth","Penta")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

DBH_roots<-lm(master3$Dry_root_wt~master3$Tree_DBH)
summary(DBH_roots)

master<-filter(master,Tree_species!="")
DBH_tree_fig<-ggplot(master,aes(x=Tree_species,y=Tree_DBH)) +
  geom_boxplot(lwd=1.5) +
  ylab(label="DBH (cm)") +
  xlab(label="Tree Species") +
  scale_x_discrete(labels=c("Goeth"="Goethalsia","Penta"="Pentaclethra"),breaks=c("Goeth","Penta")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

tree_mod<-lm(master$Tree_DBH~master$Tree_species)
anova(tree_mod) # Not significant


# Root biomass by tree type -----------------------------------------------
master3<-filter(master,Exclusion=="3")
Tree_roots_fig<-ggplot(master3,aes(x=Tree_species,y=Dry_root_wt)) +
  geom_boxplot(lwd=1.5) +
  geom_smooth(method=lm) +
  ylab(label="Dry root biomass (mg)") +
  xlab(label="Tree Species") +
  scale_x_discrete(labels=c("Goeth"="Goethalsia","Penta"="Pentaclethra"),breaks=c("Goeth","Penta")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

DBH_roots<-lm(master3$Dry_root_wt~master3$Tree_DBH)
summary(DBH_roots)



# Distribution of DBH among study trees -----------------------------------

ggplot(data=master, aes(master$Tree_DBH,fill=Tree_species)) + 
  geom_density(alpha=0.7) +
  labs(fill="Tree Species",x="Tree DBH (cm)",y="Density") +
  scale_fill_discrete(labels=c("Goeth"="Goethalsia","Penta"="Pentaclethra")) +
  theme_classic()


# Fungi to bacteria ratios: rhizosphere manipulation (no trees yet) -------

BFdat<-read.csv("Bacteria:Fungi.csv")

BF_ratio<-BFdat$Bacteria_nmol_FAME_g_soil/BFdat$Fungi_nmol_FAME_g_soil

BFdat<-cbind(BFdat,BF_ratio)
BFdat$Rhizosphere<-as.factor(BFdat$Rhizosphere)

BF_figure<-ggplot(data=BFdat,aes(x=Rhizosphere,y=BF_ratio)) +
  geom_boxplot(lwd=1.5) +
  labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold("Bacteria:Fungi Raito"))) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  ylim(0,30) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

BFmod<-lm(BFdat$BF_ratio~BFdat$Rhizosphere)
anova(BFmod)


# Microbial stress ratios --------------------------------------------------
stress_dat<-filter(data1_wide,Rhizo_Manipulation!="b")

stress_rhizo<-ggplot(data=stress_dat,aes(x=Rhizo_Manipulation,y=Stress_Ratio)) +
  geom_boxplot(lwd=1.5) +
  labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold("Microbial Stress Ratio (Saturated/Unsatured FA)"))) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))


# Microbial community by rhizosphere manipulation -------------------------

FAME_rhizo<-ggplot(data,aes(x=Rhizosphere_Manipulation,y=umol_FAME,fill=Microbial_Group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  ylab(label="FAME concentration (\u03BCmol FAME per g dry soil)") +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  xlab(label="Rhizosphere Manipulation") +
  labs(fill="Microbial Group") +
  ylim(0,42) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

FAME_rhizo_stacked_bar<-ggplot(data,aes(x=Rhizosphere_Manipulation,y=umol_FAME,fill=Microbial_Group)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(palette = "Set1") +
  ylab(label="FAME concentration (\u03BCmol FAME per g dry soil)") +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  xlab(label="Rhizosphere Manipulation") +
  labs(fill="Microbial Group") +
  #ylim(0,42) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))


# Active microbial community (d13) ----------------------------------------

d13_group_rhizo<-ggplot(data,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  ylab(label="d13 Fatty Acids") +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  xlab(label="Rhizosphere Manipulation") +
  labs(fill="Microbial Group") +
  ylim(0,42) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

