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
  annotate("text", x = 1, y = 170, label = "A", size=8, color="red") +
  annotate("text", x = 2, y = 170, label = "B", size=8, color="red") +
  annotate("text", x = 3, y = 170, label = "B", size=8, color="red") +
  ylim(-15,170) +
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
  annotate("text", x = 1, y = 170, label = "A", size=8, color="red") +
  annotate("text", x = 2, y = 170, label = "B", size=8, color="red") +
  annotate("text", x = 3, y = 170, label = "B", size=8, color="red") +
  ylim(-15,170) +
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
# Consider removing y axis label from starch graph

#### Is there a difference between leaf and starch 13CO2 release?

# New data frame, puts d13 response (Day 4 for starch, Day 9 for leaf) in one column.
smallStarchDat<-data.frame(core_ID=starchDat$core_ID, Isotope_label=starchDat$Isotope_label, Exclusion=starchDat$Exclusion, d13_CO2=starchDat$d13_d4)

smallLeafDat<-data.frame(core_ID=leafDat$core_ID, Isotope_label=leafDat$Isotope_label, Exclusion=leafDat$Exclusion, d13_CO2=leafDat$d13_d9)

comboDat<-rbind(smallStarchDat,smallLeafDat)

# Anova with interaction term
comboMod<-lm(comboDat$d13_CO2~comboDat$Exclusion*comboDat$Isotope_label)
anovaComboMod<-aov(comboMod)
anova(comboMod)
TukeyCombo<-TukeyHSD(anovaComboMod) # 13CO2 differs by exclusion but not amongst label types. The full exclosure treatment is different from all other treatments (much greater decomposition of both starch and leaf substrates when roots and mycorrhizae are absent).

# 3L diff from 1L
# 3S diff from 1L
# 1S diff from 2L
# 1S diff from 3L
# 2S diff from 1S
# 3S diff from 1S

#### Dry root weight by core type
ggplot(dat,aes(x=Exclusion,y=Dry_root_wt)) +
  geom_boxplot()

#### Enzyme activity by rhizosphere exclusion treatment

# Create composite columns for C-degrading enzyme activity (AG, BG, CBH) and nutrient acquiriing enzyme activity (NAG, PHOS, LAP)

enzyDat<-filter(dat,a_gluc!="NA")
enzyDat<-filter(enzyDat,b_gluc!="NA")
enzyDat<-filter(enzyDat,cbh!="NA")
enzyDat<-filter(enzyDat,nag!="NA")
enzyDat<-filter(enzyDat,phos!="NA")
enzyDat<-filter(enzyDat,lap!="NA")

C_enzymes<-(enzyDat$a_gluc+ enzyDat$b_gluc+ enzyDat$cbh)
nut_enzymes<-(enzyDat$nag+ enzyDat$phos+ enzyDat$lap)

enzyDat<-cbind(enzyDat,C_enzymes,nut_enzymes)

# ANOVA, carbon enzymes by rhizosphere manipulation treatment
CenzyMod<-lm(enzyDat$C_enzymes~enzyDat$Exclusion)
anovaCenzyMod<-aov(CenzyMod)
anova(CenzyMod)
TukeyEnzyC<-TukeyHSD(anovaCenzyMod) # No differences between C-acquiring enzyme activity between rhizosphere treatments

# Boxplots for C enzyme activity by rhizosphere manipulation
rhizo_enzyC<-ggplot(data=enzyDat,aes(x=Exclusion,y=C_enzymes)) +
  geom_boxplot(lwd=1.5) +
  labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("C Enzyme Activity")))) +
  #ylim(-15,170) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

# ANOVA, nutrient enzymes by rhizosphere manipulation treatment
NutEnzyMod<-lm(enzyDat$nut_enzymes~enzyDat$Exclusion)
anovaNutEnzyMod<-aov(NutEnzyMod)
anova(NutEnzyMod)
TukeyEnzyNut<-TukeyHSD(anovaNutEnzyMod) # No differences between nutrient-acquiring enzyme activity between rhizosphere treatments

# Boxplots for nutrient enzyme activity by rhizosphere manipulation
rhizo_enzyNut<-ggplot(data=enzyDat,aes(x=Exclusion,y=nut_enzymes)) +
  geom_boxplot(lwd=1.5) +
  labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("Nutrient Enzyme Activity")))) +
  #ylim(-15,170) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

#### Respiration by rhizosphere manipulation

