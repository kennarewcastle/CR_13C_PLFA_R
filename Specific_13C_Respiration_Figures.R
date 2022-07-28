#### Microbial Biomass Specific Respiration Figures
#### KER
#### Created: 07 December 2021
#### Modified: 
####        27 July 2022: Updated specific respiration calculations to reflect 13C-CO2 measured only on peak efflux day.


# Load packages -----------------------------------------------------------

library(dplyr)
library(ggplot2)
library(gridExtra)

# Read in data ------------------------------------------------------------

CO2<-read.csv("13C_CO2_Final_Calculations.csv")
CO2<-CO2[sort.list(CO2$Clean_Core_ID),]
PLFA<-read.csv("Total_13C_per_Core.csv")
PLFA<-PLFA[order(PLFA$Core_ID),]

# Check to see alignment of core IDs --------------------------------------
PLFA_Cores<-c(PLFA$Core_ID,rep(NA,times=19)) # PLFA data set is missing 19 points that are present in the CO2
cores<-data.frame("PLFA_Cores"=PLFA_Cores,"CO2_Cores"=CO2$Clean_Core_ID)

  # Eliminate the following CO2 data points that don't exist in the PLFA data set (number refers to row number in the CO2 dataset): 2, 11, 13, 19, 28, 32, 34, 37, 46, 48, 60, 61, 67, 71, 79, 83, 90, 94, 96)
  CO2<-CO2[c(1,3:10,12,14:18,20:27,29:31,33,35:36,38:45,47,49:59,62:66,68:70,72:78,80:82,84:89,91:93,95),]
  

# Combined PLFA, CO2 dataframe --------------------------------------------

data<-data.frame(PLFA,"Total_13C_CO2_ug"=CO2$Total_13C_ug)
names(data)[4]<-"Total_13C_FAME_ug"
data$Specific_Respiration_CO2_FAME<-data$Total_13C_CO2_ug/data$Total_13C_FAME_ug
#write.csv(data, file="Final_Specific_Respiration_Data.csv", row.names=FALSE)

# Create figures for specific respiration — leaf --------------------------
data<-read.csv("Final_Specific_Respiration_Data.csv")
leaf<-filter(data,Label_Type=="L")
leaf$Mesh_Type<-as.factor(leaf$Mesh_Type)

leaf_mod<-lm(Specific_Respiration_CO2_FAME~Mesh_Type,data=leaf)
anova(leaf_mod) # No significant differences

leaf_spec_resp<-ggplot(data=leaf,aes(x=Mesh_Type,y=Specific_Respiration_CO2_FAME)) +
  geom_boxplot(lwd=1) +
  geom_point(alpha=0.7,size=5) +
  labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold("Biomass Specific Respiration (13C-CO2/13C-PLFA)")),title = expression(bold("Leaf Substrate"))) +
  ylim(0,0.0050) + # This plot eliminates two crazy -R+M data points             
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA)) # This plot eliminates one crazy -R+M data point             

# Create figures for specific respiration — starch -------------------------
starch<-filter(data,Label_Type=="S")
starch$Mesh_Type<-as.factor(starch$Mesh_Type)

starch_mod<-lm(Specific_Respiration_CO2_FAME~Mesh_Type,data=starch)
anova(starch_mod) # No significant differences
TukeyHSD(aov(starch_mod))

starch_spec_resp<-ggplot(data=starch,aes(x=Mesh_Type,y=Specific_Respiration_CO2_FAME)) +
  geom_boxplot(lwd=1) +
  geom_point(alpha=0.7,size=5) +
  labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold("Biomass Specific Respiration (13C-CO2/13C-PLFA)")),title = expression(bold("Starch Substrate"))) +
  #ylim(0,0.0050) +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA)) 

grid.arrange(leaf_spec_resp,starch_spec_resp,ncol=2)
