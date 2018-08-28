# 13C CO2 Calculations for Costa Rica 2016
# April 27th, 2018
# Modified August 6th, 2018
# KER

library(dplyr)
library(ggplot2)
library(gridExtra)

dat<-read.csv("Costa Rica Master sheet annotated.csv")
dat<-filter(dat,dat$Exclusion!="NA")
dat$Exclusion<-as.factor(dat$Exclusion)

starchDat<-filter(dat,dat$Isotope_label=="S") # Only cores labeled with starch 13C
#starchDat<-filter(starchDat,starchDat$d13_d5<150)
starchDat<-filter(starchDat,starchDat$d13_d5!="NA")
starchDat<-filter(starchDat,starchDat$d13_d2!="NA")
starchDat<-filter(starchDat,starchDat$d13_d3!="NA")
starchDat<-filter(starchDat,starchDat$d13_d4!="NA")

leafDat<-filter(dat,dat$Isotope_label=="L") # Only cores labeled with leaf 13C
leafDat<-filter(dat,dat$d13_d4!="NA")
leafDat<-filter(dat,dat$d13_d6!="NA")
leafDat<-filter(dat,dat$d13_d7!="NA")
leafDat<-filter(dat,dat$d13_d9!="NA")


#### 13CO2 trends over 5 days, starch
d13_starch_dat<-c(starchDat$d13_d2,starchDat$d13_d3,starchDat$d13_d4,starchDat$d13_d5)
N<-length(starchDat[,1])
d13_starch_days<-c(rep(2,times=N),rep(3,times=N),rep(4,times=N),rep(5,times=N))
d13_peak_starch<-data.frame(day=d13_starch_days,d13_CO2=d13_starch_dat)

plot.new()
ggplot(data=d13_peak_starch,aes(x=day,y=d13_CO2)) +
  geom_point(colour="seagreen") +
  geom_smooth(colour="seagreen",size=1.5) +
  labs(x="Days",y="delta 13C") +
  theme(panel.grid.minor=element_blank(),axis.text=element_text(colour="black",size=14),axis.title=element_text(size=16,face="bold")) # Seems like day 2 is the peak for starch decomp?

d13_day_starch_mod<-lm(d13_peak_starch$d13_CO2~as.factor(d13_peak_starch$day))
anova(d13_day_starch_mod)
starch_peak<-aov(d13_day_starch_mod)
TukeyHSD(starch_peak)

#### 13CO2 trends over 9 days, leaf
d13_leaf_dat<-c(leafDat$d13_d4,leafDat$d13_d6,leafDat$d13_d7,leafDat$d13_d9)
N<-length(leafDat[,1])
d13_leaf_days<-c(rep(4,times=N),rep(6,times=N),rep(7,times=N),rep(9,times=N))
d13_leaf_days<-as.factor(d13_leaf_days)
d13_peak_leaf<-data.frame(day=d13_leaf_days,d13_CO2=d13_leaf_dat)

plot.new()
ggplot(data=d13_peak_leaf,aes(x=day,y=d13_CO2,group=day)) +
  #geom_point(colour="seagreen") +
  #geom_smooth(colour="seagreen",size=1.5) +
  geom_boxplot() +
  scale_x_discrete() +
  labs(x="Days",y="delta 13C") +
  theme(panel.grid.minor=element_blank(),axis.text=element_text(colour="black",size=14),axis.title=element_text(size=16,face="bold")) # Pretty even decomp of leaf substrate, just go with harvest day reading.

d13_day_leaf_mod<-lm(d13_peak_leaf$d13_CO2~d13_peak_leaf$day)
anova(d13_day_leaf_mod) # No statistically significant difference between days

#### 13C CO2 leaf boxplot figure
leaf_d13<-ggplot(data=leafDat,aes(x=Exclusion,y=d13_d9)) +
  geom_boxplot(lwd=1.5) +
  labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("\u03B4"^{bold("13")},"CO"[bold("2")]," (\u2030)"))),title = expression(bold("Leaf Substrate"))) +
  annotate("text", x = 1, y = 100, label = "A", size=8, color="red") +
  annotate("text", x = 2, y = 100, label = "B", size=8, color="red") +
  annotate("text", x = 3, y = 100, label = "B", size=8, color="red") +
  ylim(-20,100) +
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
starch_d13<-ggplot(data=starchDat,aes(x=Exclusion,y=d13_d2)) +
  geom_boxplot(lwd=1.5) +
  labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("\u03B4"^{bold("13")},"CO"[bold("2")]," (\u2030)"))),title = expression(bold("Starch Substrate"))) +
  annotate("text", x = 1, y = 260, label = "A", size=8, color="red") +
  annotate("text", x = 2, y = 260, label = "B", size=8, color="red") +
  annotate("text", x = 3, y = 260, label = "B", size=8, color="red") +
  ylim(-20,260) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

#### 13C CO2 starch anova and TukeyHSD
starchMod<-lm(starchDat$d13_d2~starchDat$Exclusion) # super super significant p < 0.001
anovaStarchMod<-aov(starchMod)
anova(starchMod)
TukeyHSD(anovaStarchMod) # Mesh 1 is different from the others

#### Panelled figure
grid.arrange(leaf_d13,starch_d13,nrow=1)


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
  labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("C Enzyme Activity (nmol ","g dry ","soil"^{bold("-1")}," h"^{bold("-1")},")")))) +
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
  labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("Nutrient Enzyme Activity (nmol ","g dry ","soil"^{bold("-1")}," h"^{bold("-1")},")")))) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))


#### Panelled figure
grid.arrange(rhizo_enzyC,rhizo_enzyNut,nrow=1)

#### Respiration by rhizosphere manipulation

# Need new column that lists harvest day respiration measurements for both starch (day 4) and leaf cores

#########################################################################################################
# FUNCTION: CollectResp
# Pulls the correct respiration data point corresponding with core harvest date.
# input: master data frame
# output: Data frame with core ID, exclusion type, tree type, respiration rate on harvest day
#--------------------------------------------------------------------------------------------------------

CollectResp<-function(dat) {
  dat<-filter(dat,Isotope_label!="SKIP")
  N<-length(dat[,1])
  resp<-c()
  
  for (i in 1:N) {
    if (dat[i,4]=="L") {
      resp[i]<-dat[i,35]
    }
    if (dat[i,4]=="LW") {
      resp[i]<-dat[i,35]
    }
    if (dat[i,4]=="S") {
      resp[i]<-dat[i,32]
    }
    if (dat[i,4]=="SW") {
      resp[i]<-dat[i,32]
    }
  }
  
  resp_data<-data.frame(core_ID=dat[,1],Exclusion=as.factor(dat[,2]),Tree_species=dat[,3],resp=resp)
  return(resp_data)
}

#--------------------------------------------------------------------------------------------------------

respiration_data<-CollectResp(dat=dat)

# ANOVA for respiration by rhizosphere manipulation
resp_mod<-lm(respiration_data$resp~respiration_data$Exclusion) # Marginally significant p = 0.059
anova(resp_mod)
resp_mod_aov<-aov(resp_mod)
TukeyHSD(resp_mod_aov)

# Boxplot for respiration by rhizosphere manipulation 
resp_rhizo<-ggplot(data=respiration_data,aes(x=Exclusion,y=resp)) +
  geom_boxplot(lwd=1.5) +
  labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("Respiration Rate")))) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

#### Respiration by tree type
resp_tree<-lm(respiration_data$resp~respiration_data$Tree_species) # Not significant
anova(resp_tree)

# Boxplot
resp_tree_box<-ggplot(data=respiration_data,aes(x=Tree_species,y=resp,group=Tree_species)) +
  geom_boxplot(lwd=1.5) +
  labs(x=expression(bold("Tree Species")),y=expression(bold(paste("Respiration Rate")))) +
  # scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))
