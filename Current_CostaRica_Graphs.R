#### Current Costa Rica Figures
#### KER
#### Created on: 07 December 2021
#### Last modified:
####          10 December 2021: Modified figures to include tree differences

# Load packages -----------------------------------------------------------

library(dplyr)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(viridis)
library(car)
library(agricolae)

# Color pallette for all 2 fill figures -----------------------------------

my_cols<-viridis(n=4)[2:3]

# delta13C-CO2 ------------------------------------------------------------

  ## Read in data
  dat<-read.csv("Costa Rica Master sheet annotated.csv")
  dat<-filter(dat,dat$Exclusion!="NA")
  dat$Exclusion<-as.factor(dat$Exclusion)

  ## Separate starch and leaf data
  starch_CO2<-filter(dat,dat$Isotope_label=="S") 
  leaf_CO2<-filter(dat,dat$Isotope_label=="L")

  ## Leaf d13 figure
  leaf_d13<-ggplot(data=leaf_CO2,aes(x=Exclusion,y=d13_d9,fill=Tree_species)) +
    geom_dotplot(binaxis="y",stackdir="center",alpha=0.7,position=position_dodge(0.8)) +
    geom_boxplot(lwd=1,outlier.shape=NA) +
    labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("\u03B4"^{bold("13")},"CO"[bold("2")]," (\u2030)"))),title = expression(bold("Leaf Substrate")),fill=expression(bold("Tree Species"))) +
    ylim(-20,100) +
    scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
    scale_fill_manual(values=my_cols,labels=c("G. meintha","P. macroloba")) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12,face="italic"),
        legend.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA),
        legend.position="none")

  ## 13C CO2 leaf anova and TukeyHSD
  leafMod<-lm(d13_d9~Exclusion*Tree_species,data=leaf_CO2) # super super significant p < 0.001
  Anova(leafMod,type="II") # Significant rhizosphere manipulation, tree, and tree*rhizo interaction
  TukeyHSD(aov(leafMod)) 
  tx<-with(leaf_CO2, interaction(Exclusion, Tree_species))
  amod<-aov(d13_d9 ~ tx, data=leaf_CO2)
  HSD.test(amod, "tx", group=TRUE, console=TRUE)
  
  ## Starch d13 figure
  starch_d13<-ggplot(data=starch_CO2,aes(x=Exclusion,y=d13_d2,fill=Tree_species)) +
    geom_dotplot(binaxis="y",stackdir="center",alpha=0.7,position=position_dodge(0.8)) +
    geom_boxplot(lwd=1,outlier.shape=NA) +
    labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("\u03B4"^{bold("13")},"CO"[bold("2")]," (\u2030)"))),title = expression(bold("Starch Substrate")),fill=expression(bold("Tree Species"))) +
    ylim(-20,260) +
    scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
    scale_fill_manual(values=my_cols,labels=c("G. meintha","P. macroloba")) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          axis.text.y=element_text(colour="black",size=10),
          axis.text.x=element_text(colour="black",size=14),
          axis.title=element_text(size=14,face="bold"),
          legend.text=element_text(size=12,face="italic"),
          legend.title=element_text(size=14,face="bold"),
          panel.border=element_rect(fill=NA,colour="black",size=1.5),
          panel.background=element_rect(fill=NA),
          legend.position="none")

  ## 13C CO2 starch anova and TukeyHSD
  starchMod<-lm(d13_d2~Exclusion*Tree_species,data=starch_CO2) # super super significant p < 0.001
  Anova(starchMod,type="II") # Significant rhizosphere manipulation
  TukeyHSD(aov(starchMod))
  tx<-with(starch_CO2, interaction(Exclusion, Tree_species))
  amod<-aov(d13_d2 ~ tx, data=starch_CO2)
  HSD.test(amod, "tx", group=TRUE, console=TRUE)

  ## Panelled figure
  grid.arrange(leaf_d13,starch_d13,nrow=1)


# delta13C-PLFA -----------------------------------------------------------
  
  ## Read in data
  data<-read.csv("FINAL_PLFA_with_metadata.csv")
  data<-data[,2:15] # Gets rid of ID column
  data<-filter(data,Microbial_Group!="standard")
  data$Rhizosphere_Manipulation<-as.factor(data$Rhizosphere_Manipulation)
  data$umol_FAME<-data$nmol_FAME_per_g_soil/1000
  
  ## Simplify the microbial groupings
  data_trim<-data.frame("Sample_Name"=data[,1],"Tree"=data[,2],"Rhizosphere_Manipulation"=data[,3],"Substrate_Type"=data[,4],"Microbial_Group"=data[,14],"d13"=data[,9])
  
  data_wide<-dcast(data_trim,Sample_Name+Tree+Rhizosphere_Manipulation+Substrate_Type ~ Microbial_Group, value.var="d13", fun.aggregate=mean, na.rm=TRUE)
  
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))
  data_wide[is.nan(data_wide)] <- NA
  
      ### Take mean of narrow groups to yield FUNGI, BACTERIA, ALL_MICROBE groups
      for (i in 1:nrow(data_wide)) {
        data_wide$BACTERIA[i]<-mean(c(data_wide$bacteria[i],data_wide$`gram+bacteria`[i],data_wide$`gram-bacteria`[i]),na.rm=TRUE)
        data_wide$FUNGI[i]<-mean(c(data_wide$AMF[i],data_wide$fungi[i]),na.rm=TRUE)
        data_wide$ALL_MICROBES[i]<-mean(c(data_wide$AMF[i],data_wide$bacteria[i],data_wide$fungi[i],data_wide$general[i],data_wide$`gram-bacteria`[i],data_wide$`gram+bacteria`[i],data_wide$protozoa[i]),na.rm=TRUE)
      }

  data_wide<-data_wide[,c(1,2,3,4,5,12,13,14)] # Trim down to broader groups
  
  data<-melt(data_wide,value.name="d13",id.vars=c("Sample_Name","Rhizosphere_Manipulation","Tree","Substrate_Type"))
  names(data)[5]<-"Microbial_Group"
  
  ## Add in value to create shading for rhizosphere treatments
  for (i in 1:nrow(data)) {
    if (data$Rhizosphere_Manipulation[i]=="1" | data$Rhizosphere_Manipulation[i]=="3") {
      data$stripe[i]<-"#33333333"
    } else {data$stripe[i]<-"#00000000"}
  }

  ## Separate leaf and starch datasets
  leaf_PLFA<-filter(data,Substrate_Type=="L")
  starch_PLFA<-filter(data,Substrate_Type=="S")  

  ## Leaf PLFA figure
  tree_labs<-c("G. meintha", "P. macroloba")
  names(tree_labs)<-c("G", "P")
  
  leaf_d13_PLFA<-ggplot(leaf_PLFA,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
    geom_rect(xmin=1.5,xmax=2.5,ymin=-Inf,ymax=Inf,fill="gray70") +
    geom_hline(yintercept=-25,lwd=1.5,linetype="dashed") +
    geom_boxplot(lwd=0.75) +
    ylab(label=expression(bold(paste("\u03B4"^{bold("13")}, "C PLFA"," (\u2030)")))) +
    scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
    scale_fill_viridis(discrete=TRUE,labels=c("AMF","Bacteria","Fungi","All Microbes")) +
    xlab(label="Rhizosphere Manipulation") +
    ggtitle("Leaf Substrate") +
    labs(fill="Microbial Group") +
    ylim(-50,150) +
    facet_wrap(~Tree,labeller=labeller(Tree=tree_labs)) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          axis.text.y=element_text(colour="black",size=10),
          axis.text.x=element_text(colour="black",size=14,angle=90),
          axis.title=element_text(size=14,face="bold"),
          plot.title=element_text(size=14,face="bold"),
          panel.border=element_rect(fill=NA,colour="black",size=1.5),
          panel.background=element_rect(fill=NA),
          strip.text.x=element_text(size=14,face="bold.italic"),
          legend.position="none")
  
  ## Leaf PLFA ANOVA and Tukey's HSD
  leaf_PLFA_mod<-lm(data=leaf_PLFA,d13~Rhizosphere_Manipulation*Microbial_Group*Tree)
  Anova(leaf_PLFA_mod,type="II") # Rhizosphere manipulation significant (p < 0.001), rhizosphere x tree significant
  TukeyHSD(aov(leaf_PLFA_mod))
  tx<-with(leaf_PLFA, interaction(Rhizosphere_Manipulation,Microbial_Group,Tree))
  amod<-aov(d13 ~ tx, data=leaf_PLFA)
  HSD.test(amod, "tx", group=TRUE, console=TRUE)
  
  ## Starch PLFA figure
  tree_labs<-c("G. meintha", "P. macroloba")
  names(tree_labs)<-c("G", "P")
  
  starch_d13_PLFA<-ggplot(data=starch_PLFA,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
    geom_rect(xmin=1.5,xmax=2.5,ymin=-Inf,ymax=Inf,fill="gray70") +
    geom_hline(yintercept=-25,lwd=1.5,linetype="dashed") +
    geom_boxplot(lwd=0.75) +
    scale_fill_viridis(discrete=TRUE,labels=c("AMF","Bacteria","Fungi","All Microbes")) +
    ylab(label=expression(bold(paste("\u03B4"^{bold("13")}, "C PLFA"," (\u2030)")))) +
    scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
    xlab(label="Rhizosphere Manipulation") +
    geom_hline(yintercept=-25,lwd=1.5,linetype="dashed") +
    ggtitle("Starch Substrate") +
    labs(fill="Microbial Group") +
    ylim(-50,150) + # Excludes 1 outlier
    facet_wrap(~Tree,labeller=labeller(Tree=tree_labs)) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          axis.text.y=element_text(colour="black",size=10),
          axis.text.x=element_text(colour="black",size=14,angle=90),
          axis.title=element_text(size=14,face="bold"),
          plot.title=element_text(size=14,face="bold"),
          panel.border=element_rect(fill=NA,colour="black",size=1.5),
          panel.background=element_rect(fill=NA),
          strip.text.x=element_text(size=14,face="bold.italic"),
          legend.position="none")
  
  ## Starch PLFA ANOVA and Tukey's HSD
  starch_PLFA_mod<-lm(data=starch_PLFA,d13~Rhizosphere_Manipulation*Microbial_Group*Tree)
  Anova(starch_PLFA_mod,type="II") 
  TukeyHSD(aov(starch_PLFA_mod))
  tx<-with(starch_PLFA, interaction(Rhizosphere_Manipulation,Microbial_Group,Tree))
  amod<-aov(d13 ~ tx, data=starch_PLFA)
  HSD.test(amod, "tx", group=TRUE, console=TRUE)
  
  ## Panneled figure
  grid.arrange(leaf_d13_PLFA,starch_d13_PLFA,ncol=2)
  

# 13C Specific Respiration ------------------------------------------------

  ## Read in data
  data<-read.csv("Final_Specific_Respiration_Data.csv")
  data$Mesh_Type<-as.factor(data$Mesh_Type)
  leaf<-filter(data,Label_Type=="L")
  starch<-filter(data,Label_Type=="S")

  ## Leaf specific respiration fig
  leaf_spec_resp<-ggplot(data=leaf,aes(x=Mesh_Type,y=Specific_Respiration_CO2_FAME,fill=Tree_Species)) +
    geom_dotplot(binaxis="y",stackdir="center",alpha=0.7,position=position_dodge(0.8)) +
    geom_boxplot(lwd=1,outlier.shape=NA) +
    labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("Biomass Specific Respiration ","("^{bold("13")},"C-","CO"[bold("2")]," / "^{bold("13")},"C-PLFA)"))),title = expression(bold("Leaf Substrate")),fill=expression(bold("Tree Species"))) +
    ylim(0,0.0050) + # This plot eliminates two crazy -R+M data points             
    scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
    scale_fill_manual(values=my_cols,labels=c("G. meintha","P. macroloba")) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          axis.text.y=element_text(colour="black",size=10),
          axis.text.x=element_text(colour="black",size=14),
          axis.title=element_text(size=14,face="bold"),
          legend.text=element_text(size=12,face="italic"),
          legend.title=element_text(size=14,face="bold"),
          panel.border=element_rect(fill=NA,colour="black",size=1.5),
          panel.background=element_rect(fill=NA),
          legend.position="none")
  
  ## Leaf specific respiration ANOVA and TukeyHSD
  leaf_mod<-lm(Specific_Respiration_CO2_FAME~Mesh_Type*Tree_Species,data=leaf)
  Anova(leaf_mod,type="II")
  TukeyHSD(aov(starch_PLFA_mod))
  tx<-with(leaf, interaction(Mesh_Type,Tree_Species))
  amod<-aov(Specific_Respiration_CO2_FAME ~ tx, data=leaf)
  HSD.test(amod, "tx", group=TRUE, console=TRUE)
  
  ## Starch specific respiration fig
  starch_spec_resp<-ggplot(data=starch,aes(x=Mesh_Type,y=Specific_Respiration_CO2_FAME,fill=Tree_Species)) +
    geom_dotplot(binaxis="y",stackdir="center",alpha=0.7,position=position_dodge(0.8)) +
    geom_boxplot(lwd=1,outlier.shape=NA) +
    labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("Biomass Specific Respiration ","("^{bold("13")},"C-","CO"[bold("2")]," / "^{bold("13")},"C-PLFA)"))),title = expression(bold("Starch Substrate")),fill=expression(bold("Tree Species"))) +
    scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
    scale_fill_manual(values=my_cols,labels=c("G. meintha","P. macroloba")) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          axis.text.y=element_text(colour="black",size=10),
          axis.text.x=element_text(colour="black",size=14),
          axis.title=element_text(size=14,face="bold"),
          legend.text=element_text(size=12,face="italic"),
          legend.title=element_text(size=14,face="bold"),
          panel.border=element_rect(fill=NA,colour="black",size=1.5),
          panel.background=element_rect(fill=NA),
          legend.position="none")
  
  ## Starch specific respiration ANOVA and TukeyHSD
  starch_mod<-lm(Specific_Respiration_CO2_FAME~Mesh_Type*Tree_Species,data=starch)
  Anova(starch_mod,type="II") # No significant differences
  TukeyHSD(aov(starch_mod))
  tx<-with(starch, interaction(Mesh_Type,Tree_Species))
  amod<-aov(Specific_Respiration_CO2_FAME ~ tx, data=starch)
  HSD.test(amod, "tx", group=TRUE, console=TRUE)

# Respiration (non-labelled) ----------------------------------------------

  ## Read in data
  data<-read.csv("Costa Rica Master sheet annotated.csv")
  data<-filter(data,Exclusion!="NA") # Several mislabeled cores for which we can't determine rhizosphere treatment; exclude these
  data$Exclusion<-as.factor(data$Exclusion)
  
  ## Create Mean_Respiration variable that averages respiration rates across the 9 sample days (different cores measured respiration on different days, averaging would even out some of the differences in respiration due to the sampling schedule)
  for (i in 1:nrow(data)) {
    data$Mean_Respiration[i]<-mean(c(data$Efflux_d1[i],data$Efflux_d2[i],data$Efflux_d3[i],data$Efflux_d4[i],data$Efflux_d5[i],data$Efflux_d6[i],data$Efflux_d7[i],data$Efflux_d9[i]),na.rm=TRUE)
  }
  
  data<-filter(data,Mean_Respiration!="NaN") # Remove NaNs produced from entries that have no respiration data
  #write.csv(data,file="Resp_Data_for_Spec_Resp.csv",row.names=FALSE)
  
  ## Respiration figure
  resp_fig<-ggplot(data=data,aes(x=Exclusion,y=Mean_Respiration,fill=Tree_species)) +
    geom_dotplot(binaxis="y",stackdir="center",alpha=0.7,position=position_dodge(0.8)) +
    geom_boxplot(lwd=1,outlier.shape=NA) +
    labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("Respiration Rate (",mu,"mol C",O[bold("2")]," ",m^{bold("-2")},s^{bold("-1")},")"))),title = expression(bold("Soil Respiration")),fill="Tree Species") +
    scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
    scale_fill_manual(values=my_cols,labels=c("G. meintha","P. macroloba")) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          axis.text.y=element_text(colour="black",size=10),
          axis.text.x=element_text(colour="black",size=14),
          axis.title=element_text(size=14,face="bold"),
          legend.title=element_text(size=14,face="bold"),
          legend.text=element_text(size=12,face="italic"),
          panel.border=element_rect(fill=NA,colour="black",size=1.5),
          panel.background=element_rect(fill=NA),
          legend.position="none")
  
  ## Soil respiration ANOVA and Tukey HSD
  resp_mod<-lm(data=data,Mean_Respiration~Exclusion*Tree_species)
  Anova(resp_mod,type="II") # Exclusion significant
  TukeyHSD(aov(resp_mod))
  tx<-with(data, interaction(Exclusion,Tree_species))
  amod<-aov(Mean_Respiration ~ tx, data=data)
  HSD.test(amod, "tx", group=TRUE, console=TRUE)

# Unlabeled PLFA ----------------------------------------------------------

  ## Read in data
  data<-read.csv("FINAL_PLFA_with_metadata.csv")
  data<-data[,2:15] # Gets rid of ID column
  data<-filter(data,Microbial_Group!="standard")
  data$Rhizosphere_Manipulation<-as.factor(data$Rhizosphere_Manipulation)
  data$umol_FAME<-data$nmol_FAME_per_g_soil/1000
  
  ## Simplify the microbial groupings
  data_trim<-data.frame("Sample_Name"=data[,1],"Tree"=data[,2],"Rhizosphere_Manipulation"=data[,3],"Microbial_Group"=data[,14],"umol_FAME"=data[,15])
  
  data_wide<-dcast(data_trim,Sample_Name+Tree+Rhizosphere_Manipulation ~ Microbial_Group, value.var="umol_FAME", fun.aggregate=sum, na.rm=TRUE)
  
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))
  data_wide[is.nan(data_wide)] <- NA
  
  ### Add together mass of FAME for narrow groups to yield FUNGI, BACTERIA, ALL_MICROBE groups
  for (i in 1:nrow(data_wide)) {
    data_wide$BACTERIA[i]<-sum(c(data_wide$bacteria[i],data_wide$`gram+bacteria`[i],data_wide$`gram-bacteria`[i]))
    data_wide$FUNGI[i]<-sum(c(data_wide$AMF[i],data_wide$fungi[i]))
    data_wide$ALL_MICROBES[i]<-sum(c(data_wide$AMF[i],data_wide$bacteria[i],data_wide$fungi[i],data_wide$general[i],data_wide$`gram-bacteria`[i],data_wide$`gram+bacteria`[i],data_wide$protozoa[i]))
  }
  
  data_wide<-data_wide[,c(1,2,3,4,11,12,13)] # Trim down to broader groups
  
  data<-melt(data_wide,value.name="umol_FAME",id.vars=c("Sample_Name","Tree","Rhizosphere_Manipulation"))
  names(data)[4]<-"Microbial_Group"
  data$ln_FAME<-log(data$umol_FAME)
  #write.csv(data,file="Unlab_PLFA_for_Spec_Resp.csv",row.names=FALSE)
  
  ## Un-labelled PLFA figure
  tree_labs<-c("G. meintha", "P. macroloba")
  names(tree_labs)<-c("G", "P")
  
  reg_PLFA_fig<-ggplot(data=data,aes(x=Rhizosphere_Manipulation,y=ln_FAME,fill=Microbial_Group)) +
    geom_rect(xmin=1.5,xmax=2.5,ymin=-Inf,ymax=Inf,fill="gray70") +
    geom_boxplot(lwd=0.75) +
    scale_fill_viridis(discrete=TRUE,labels=c("AMF","Bacteria","Fungi","All Microbes")) +
    ylab(expression(bold(paste("FAME Concentration (ln(",mu,"mol FAME ",g^{bold("-1")}," dry soil))")))) +
    scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
    xlab(label="Rhizosphere Manipulation") +
    ggtitle("Background PLFA") +
    labs(fill="Microbial Group") +
    facet_wrap(~Tree,labeller=labeller(Tree=tree_labs)) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          axis.text.y=element_text(colour="black",size=10),
          axis.text.x=element_text(colour="black",size=14,angle=90),
          axis.title=element_text(size=14,face="bold"),
          plot.title=element_text(size=14,face="bold"),
          panel.border=element_rect(fill=NA,colour="black",size=1.5),
          panel.background=element_rect(fill=NA),
          strip.text.x=element_text(size=14,face="bold.italic"),
          legend.position="none")
  
  ## ANOVA for general PLFA by rhizosphere manipulation and tree type
  PLFA_mod<-lm(data=data,umol_FAME~Rhizosphere_Manipulation*Microbial_Group*Tree)
  Anova(PLFA_mod,type="II") # Rhizospehre manipulation not significant, microbial group significant, tree significant, rhizosphere*microbial group not significant
  

# Specific respiration (unlabelled) ---------------------------------------

  ## Read in data
  spec_PLFA<-read.csv("Unlab_PLFA_for_Spec_Resp.csv") # 77 observations
  spec_PLFA<-filter(spec_PLFA,Microbial_Group=="ALL_MICROBES")
  
  spec_CO2<-read.csv("Resp_Data_for_Spec_Resp.csv") # 138 observations
  
  ## Eliminate CO2 cores that are not present in the PLFA dataset
  PLFA_Cores<-c(spec_PLFA$Sample_Name,rep(NA,times=61)) # PLFA data set is missing 61 observations
  cores<-data.frame("PLFA_Cores"=PLFA_Cores,"CO2_Cores"=spec_CO2$core_ID)
  
    ### Eliminate the following CO2 data points that doint exist in the PLFA data set (number refers to row number in the CO2 dataset): 2, 7, 8, 9, 14, 16, 17, 18, 19, 22, 23, 24, 28, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 49, 50, 51, 52, 56, 58, 61, 67, 68, 69, 73, 75, 85, 86, 87, 90, 91, 92, 93, 94, 95, 96, 97, 103, 106, 107, 108, 110, 118, 119, 120, 121, 125, 132, 136, 138 
    CO2<-spec_CO2[c(1,3:6,10:13,15,20:21,25:27,29:33,46:48,53:55,57,59:60,62:66,70:72,74,76:84,88:89,98:102,104:105,109,111:117,122:124,126:131,133:135,137),]
  
  ## Combine CO2 and PLFA data frames
  unlab_data<-data.frame(spec_PLFA,"Mean_Respiration"=CO2$Mean_Respiration)
  unlab_data$Specific_Respiration_CO2_FAME<-unlab_data$Mean_Respiration/unlab_data$umol_FAME
  unlab_data$Rhizosphere_Manipulation<-as.factor(unlab_data$Rhizosphere_Manipulation)
  
  ## Unlabelled specific respiration figure
  unlab_spec_resp<-ggplot(data=unlab_data,aes(x=Rhizosphere_Manipulation,y=Specific_Respiration_CO2_FAME,fill=Tree)) +
    geom_dotplot(binaxis="y",stackdir="center",alpha=0.7,position=position_dodge(0.8)) +
    geom_boxplot(lwd=1,outlier.shape=NA) +
    labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("Specific Respiration ","(",mu,"mol ","CO"[bold("2")]," ",m^{bold("-2")}," ",s^{bold("-1")}," / ",mu,"mol FAME ",g^{bold("-1")}," dry soil"))),title = expression(bold("Unlabelled Carbon")),fill=expression(bold("Tree Species"))) +
    scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
    scale_fill_manual(values=my_cols,labels=c("G. meintha","P. macroloba")) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          axis.text.y=element_text(colour="black",size=10),
          axis.text.x=element_text(colour="black",size=14),
          axis.title.x=element_text(size=14,face="bold"),
          axis.title.y=element_text(size=12,face="bold"),
          legend.text=element_text(size=12,face="italic"),
          legend.title=element_text(size=14,face="bold"),
          panel.border=element_rect(fill=NA,colour="black",size=1.5),
          panel.background=element_rect(fill=NA),
          legend.position="none")
  
  ## Unlabeled specific respiration ANOVA and TukeyHSD
  unlab_mod<-lm(Specific_Respiration_CO2_FAME~Rhizosphere_Manipulation*Tree,data=unlab_data)
  Anova(unlab_mod,type="II") # No significant differences
  TukeyHSD(aov(unlab_mod))
  tx<-with(unlab_data, interaction(Rhizosphere_Manipulation,Tree))
  amod<-aov(Specific_Respiration_CO2_FAME ~ tx, data=unlab_data)
  HSD.test(amod, "tx", group=TRUE, console=TRUE)
 
# Pannels of figures with labelled leaf and starch, unlabelled C ----------

  grid.arrange(leaf_d13,starch_d13,resp_fig,nrow=1) # CO2
  grid.arrange(leaf_d13_PLFA,starch_d13_PLFA,reg_PLFA_fig,nrow=1) # PLFA
  grid.arrange(leaf_spec_resp,starch_spec_resp,unlab_spec_resp,nrow=1) # Specific respiration
  
  

  
  
  