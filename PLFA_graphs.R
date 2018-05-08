# Constructing graphs for PLFA data
# 8 May 2018
# KER


# Preliminaries
library(ggplot2)
library(dplyr)
library(gridExtra)

data<-read.csv("FINAL_PLFA_with_metadata.csv")
names(data)

# Remove standards from data set
data<-filter(data,Microbial_Group!="standard")

# Bar graph of microbial groups by tree type
ggplot(data,aes(x=Tree_Type,y=nmol_FAME_per_g_soil,fill=Microbial_Group)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(palette = "Set2") # Looks like no differences

# Bar graph of microbial groups by rhizosphere manipulation
ggplot(data,aes(x=Rhizosphere_Manipulation,y=nmol_FAME_per_g_soil,fill=Microbial_Group)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(palette = "Set2") # All groups increase with mesh size?

# Bar graph of d13 microbial groups by tree type
ggplot(data,aes(x=Tree_Type,y=d13,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_brewer(palette = "Set2") +
  ylab(label="d13 FAME compounds")

# Bar graph of d13 microbial groups by rhizosphere treatment
ggplot(data,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_brewer(palette = "Set2") +
  ylab(label="d13 FAME compounds")

# Splitting bar graph of rhizosphere manipulation into starch and leaf
starch_dat<-filter(data,Substrate_Type=="S") # filter data
leaf_dat<-filter(data,Substrate_Type=="L")

starch_rhizo<-ggplot(starch_dat,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_brewer(palette = "Set1") +
  ggtitle("Starch Substrate") +
  ylab(label="d13 FAME compounds")

leaf_rhizo<-ggplot(leaf_dat,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_brewer(palette = "Set1") +
  ggtitle("Leaf Substrate") +
  ylab(label="d13 FAME compounds")
  

grid.arrange(starch_rhizo,leaf_rhizo,nrow=1) # THE COOLEST GRAPH OF THEM ALL!!!

