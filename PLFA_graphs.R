# Constructing graphs for PLFA data
# 7 May 2018
# Edited: 27 September 2018
# KER


# Preliminaries
library(ggplot2)
library(dplyr)
library(gridExtra)

data<-read.csv("FINAL_PLFA_with_metadata.csv")
data<-filter(data,Microbial_Group!="standard") # Remove standards from data set
names(data)
data$Rhizosphere_Manipulation<-as.factor(data$Rhizosphere_Manipulation)
umol_FAME<-data$nmol_FAME_per_g_soil/1000
data<-cbind(data,umol_FAME)

# Bar graph of microbial groups by tree type
ggplot(data,aes(x=Tree_Type,y=umol_FAME,fill=Microbial_Group)) +
  geom_bar(stat="identity") +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() # Looks like no differences

ggplot(data,aes(x=Tree_Type,y=umol_FAME,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() # Looks like no differences

# Box plot of microbial groups by tree type

# Bar graph of microbial groups by rhizosphere manipulation
# Two panneled figure, one for each tree type.

data_goeth<-filter(data,Tree_Type=="G")
FAME_rhizo_g<-ggplot(data_goeth,aes(x=Rhizosphere_Manipulation,y=umol_FAME,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
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
        panel.background=element_rect(fill=NA),
        legend.position="none")


data_pent<-filter(data,Tree_Type=="P")
FAME_rhizo_p<-ggplot(data_pent,aes(x=Rhizosphere_Manipulation,y=umol_FAME,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
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

library(gridExtra)
grid.arrange(FAME_rhizo_p,FAME_rhizo_g,nrow=1)

#ggsave(filename="FAME_rhizosphere.jpg",plot=FAME_rhizo)

FAME_rhizo_mod<-lm(data$umol_FAME~data$Rhizosphere_Manipulation*data$Microbial_Group) # super super significant p < 0.001
anovaFAME_rhizo<-aov(FAME_rhizo_mod)
anova(FAME_rhizo_mod)
TukeyHSD(anovaFAME_rhizo) # None of the microbial groups change by rhizosphere treatment


# Bar graph of d13 microbial groups by tree type
ggplot(data,aes(x=Tree_Type,y=d13,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_brewer(palette = "Set1") +
  ylab(label="d13 FAME compounds") +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

# Box plot of d13 microbial groups by tree type
ggplot(data,aes(x=Tree_Type,y=d13,fill=Microbial_Group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  ylab(label="d13 FAME compounds")

# Bar graph of d13 microbial groups by rhizosphere treatment

starch13<-filter(data,Substrate_Type=="S")
leaf13<-filter(data,Substrate_Type=="L")

starch13_plot<-ggplot(starch13,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_brewer(palette = "Set1") +
  ylab(label="d13 FAME compounds") +
  ggtitle("Starch Label") +
  ylim(-100,1000) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  xlab(label="Rhizosphere Manipulation") +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA),
        legend.position="none")

leaf13_plot<-ggplot(leaf13,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
  scale_fill_brewer(palette = "Set1") +
  ylab(label="d13 FAME compounds") +
  ggtitle("Leaf Label") +
  ylim(-100,1000) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  xlab(label="Rhizosphere Manipulation") +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        plot.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA),
        legend.position="none")

grid.arrange(starch13_plot,leaf13_plot,nrow=1)


# Box plot of d13 microbial groups by rhizosphere treatment
ggplot(data,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  ylim(-50,500) +
  ylab(label="d13 FAME compounds") +
  theme_classic()

range(data$d13)
# Splitting bar graph of rhizosphere manipulation into starch and leaf
starch_dat<-filter(data,Substrate_Type=="S") # filter data
leaf_dat<-filter(data,Substrate_Type=="L")
starch_dat$Rhizosphere_Manipulation<-as.factor(starch_dat$Rhizosphere_Manipulation)
leaf_dat$Rhizosphere_Manipulation<-as.factor(leaf_dat$Rhizosphere_Manipulation)

starch_rhizo<-ggplot(starch_dat,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
  ylim(-50,1000) +
  scale_fill_brewer(palette = "Set1") +
  ggtitle("Starch Substrate") +
  ylab(label="d13 FAME Compounds") +
  xlab(label="Rhizosphere Manipulation") +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA),
        legend.position="none")
  

leaf_rhizo<-ggplot(leaf_dat,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_bar(stat="identity",position="dodge") +
  ggtitle("Leaf Substrate") +
  ylim(-50,1000) +
  xlab(label="Rhizosphere Manipulation") +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  scale_fill_brewer(palette = "Set1") +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA),
        axis.title.y=element_blank())
  

grid.arrange(starch_rhizo,leaf_rhizo,nrow=1) # THE COOLEST GRAPH OF THEM ALL!!!

# Splitting leaf and starch, box plot
starch_rhizo_box<-ggplot(starch_dat,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set1") +
  ggtitle("Starch Substrate") +
  ylab(label="d13 FAME Compounds") +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  xlab(label="Rhizosphere Manipulation") +
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
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  xlab(label="Rhizosphere Manipulation") +
  scale_fill_discrete(name="Microbial\nGroup",breaks=unique(leaf_dat$Microbial_Group),
                      labels=c("Bacteria", "Fungi", "Gram - Bacteria", "Gram + Bacteria", "Protozoa")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

grid.arrange(starch_rhizo_box,leaf_rhizo_box,nrow=1) # THE COOLEST GRAPH OF THEM ALL!!!

# ANOVA on differences between d13 signal in microbial groups x rhizosphere manipulation--LEAF
leaf_d13<-lm(leaf_dat$d13~leaf_dat$Microbial_Group*leaf_dat$Rhizosphere_Manipulation)
anova(leaf_d13)
anovaLeafd13<-aov(leaf_d13)
TukeyHSD(anovaLeafd13) # Bacteria is unique, fungi is usually unique, 3 is different from 1 & 2

# ANOVA on differences between d13 signal in microbial groups x rhizosphere manipulation--STARCH
starch_d13<-lm(starch_dat$d13~starch_dat$Microbial_Group*starch_dat$Rhizosphere_Manipulation)
anova(starch_d13)
anovaStarchd13<-aov(starch_d13)
TukeyHSD(anovaStarchd13) # Bacteria is unique, fungi is usually unique, 3 is different from 1 & 2

# AMF abundance by rhizosphere treatment
AMF_dat<-filter(data,Microbial_Group=="AMF")

plot.new()
ggplot(data=AMF_dat,aes(x=Rhizosphere_Manipulation,y=umol_FAME,fill=Tree_Type)) +
  geom_boxplot() +
  labs(x="Rhizosphere Manipulation",y="umol FAME compound/ g dry soil") +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

# AMF abundance by tree type
ggplot(data=AMF_dat,aes(x=Tree_Type,y=umol_FAME,fill=Tree_Type)) +
  geom_boxplot() +
  labs(x="Tree Species",y="umol FAME compound/ g dry soil") +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))


# ANOVA on differences between AMF abundance between rhizosphere treatments
rhizo_AMF<-lm(AMF_dat$umol_FAME~AMF_dat$Rhizosphere_Manipulation)
anova(rhizo_AMF) # p = 0.319
anovaRhizo_AMF<-aov(rhizo_AMF)
TukeyHSD(anovaRhizo_AMF) # Bacteria is unique, fungi is usually unique, 3 is different from 1 & 2

# ANOVA on differences between AMF abundance between tree types
rhizo_tree_AMF<-lm(AMF_dat$umol_FAME~AMF_dat$Rhizosphere_Manipulation*AMF_dat$Tree_Type) 
anova(rhizo_tree_AMF) # Significant difference between AMF abundance between tree types (p = 0.0236), no interaction between tree type and rhizosphere manipulation


