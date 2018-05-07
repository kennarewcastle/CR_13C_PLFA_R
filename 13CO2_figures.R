# 13C CO2 Calculations for Costa Rica 2016
# April 27th, 2018
# KER

library(dplyr)
library(ggplot2)
library(gridExtra)

dat<-read.csv("Costa Rica Master sheet annotated.csv")
dat<-filter(dat,dat$Exclusion!="NA")
dat$Exclusion<-as.factor(dat$Exclusion)

leafDat<-filter(dat,dat$Isotope_label=="L") # Only cores labeled with leaf 13C
starchDat<-filter(dat,dat$Isotope_label=="S") # Only cores labeled with starch 13C
#starchDat<-filter(starchDat,starchDat$d13_d5<150)
starchDat<-filter(starchDat,starchDat$d13_d4!="NA")


#### 13C CO2 leaf boxplot figure
leaf_d13<-ggplot(data=leafDat,aes(x=Exclusion,y=d13_d9)) +
  geom_boxplot(lwd=1.5) +
  labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("\u03B4"^{bold("13")},"CO"[bold("2")]," (\u2030)"))),title = expression(bold("Leaf Substrate"))) +
  annotate("text", x = 1, y = 100, label = "A", size=8, color="red") +
  annotate("text", x = 2, y = 100, label = "B", size=8, color="red") +
  annotate("text", x = 3, y = 100, label = "B", size=8, color="red") +
  ylim(-15,100) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

#### 13C CO2 leaf anova and TukeyHSD
leafMod<-lm(leafDat$d13_d9~leafDat$Exclusion) # super super significant p < 0.001
anovaLeafMod<-aov(leafMod)
anova(leafMod)
TukeyHSD(anovaLeafMod) # Mesh 1 is different from the others

#### 13C CO2 starch boxplot figure
starch_d13<-ggplot(data=starchDat,aes(x=Exclusion,y=d13_d4)) +
  geom_boxplot(lwd=1.5) +
  labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("\u03B4"^{bold("13")},"CO"[bold("2")]," (\u2030)"))),title = expression(bold("Starch Substrate"))) +
  annotate("text", x = 1, y = 175, label = "A", size=8, color="red") +
  annotate("text", x = 2, y = 175, label = "B", size=8, color="red") +
  annotate("text", x = 3, y = 175, label = "B", size=8, color="red") +
  ylim(-25,175) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

#### 13C CO2 starch anova and TukeyHSD
starchMod<-lm(starchDat$d13_d4~starchDat$Exclusion) # super super significant p = 0.004
anovaStarchMod<-aov(starchMod)
anova(starchMod)
TukeyHSD(anovaStarchMod) # Mesh 1 is different from the others

#### Panelled figure
grid.arrange(leaf_d13,starch_d13,nrow=1)

