# Constructing graphs for PLFA data
# 8 May 2018
# KER


# Preliminaries
library(ggplot2)
library(dplyr)
library(gridExtra)

data<-read.csv("FINAL_PLFA_with_metadata.csv")
data<-filter(data,Microbial_Group!="standard") # Remove standards from data set
names(data)
data$Rhizosphere_Manipulation<-as.factor(data$Rhizosphere_Manipulation)

# Bar graph of microbial groups by tree type
ggplot(data,aes(x=Tree_Type,y=nmol_FAME_per_g_soil,fill=Microbial_Group)) +
  geom_bar() +
  scale_fill_brewer(palette = "Set2") # Looks like no differences

# Box plot of microbial groups by tree type

# Bar graph of microbial groups by rhizosphere manipulation
ggplot(data,aes(x=Rhizosphere_Manipulation,y=nmol_FAME_per_g_soil,fill=Microbial_Group)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(palette = "Set2") # All groups increase with mesh size?

# Bar graph of d13 microbial groups by tree type
ggplot(data,aes(x=Tree_Type,y=d13,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_brewer(palette = "Set1") +
  ylab(label="d13 FAME compounds")

# Box plot of d13 microbial groups by tree type
ggplot(data,aes(x=Tree_Type,y=d13,fill=Microbial_Group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  ylab(label="d13 FAME compounds")

# Bar graph of d13 microbial groups by rhizosphere treatment
ggplot(data,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_brewer(palette = "Set2") +
  ylab(label="d13 FAME compounds")

# Box plot of d13 microbial groups by rhizosphere treatment
ggplot(data,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  ylab(label="d13 FAME compounds")

# Splitting bar graph of rhizosphere manipulation into starch and leaf
starch_dat<-filter(data,Substrate_Type=="S") # filter data
leaf_dat<-filter(data,Substrate_Type=="L")
starch_dat$Rhizosphere_Manipulation<-as.factor(starch_dat$Rhizosphere_Manipulation)
leaf_dat$Rhizosphere_Manipulation<-as.factor(leaf_dat$Rhizosphere_Manipulation)

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

# Splitting leaf and starch, box plot
starch_rhizo_box<-ggplot(starch_dat,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  ggtitle("Starch Substrate") +
  ylab(label="d13 FAME compounds") +
  ylim(-100,500) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

leaf_rhizo_box<-ggplot(leaf_dat,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  ggtitle("Leaf Substrate") +
  ylab(label="d13 FAME compounds") +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

grid.arrange(starch_rhizo_box,leaf_rhizo_box,nrow=1) # THE COOLEST GRAPH OF THEM ALL!!!

# ANOVA on differences between d13 signal in microbial groups x rhizosphere manipulation--LEAF
leaf_d13<-lm(leaf_dat$d13~leaf_dat$Microbial_Group+leaf_dat$Rhizosphere_Manipulation)
anova(leaf_d13)
anovaLeafd13<-aov(leaf_d13)
TukeyHSD(anovaLeafd13) # Bacteria is unique, fungi is usually unique, 3 is different from 1 & 2

# ANOVA on differences between d13 signal in microbial groups x rhizosphere manipulation--STARCH
starch_d13<-lm(starch_dat$d13~starch_dat$Microbial_Group+starch_dat$Rhizosphere_Manipulation)
anova(starch_d13)
anovaStarchd13<-aov(starch_d13)
TukeyHSD(anovaStarchd13) # Bacteria is unique, fungi is usually unique, 3 is different from 1 & 2
