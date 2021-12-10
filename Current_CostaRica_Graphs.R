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
    geom_boxplot(lwd=1.5,outlier.shape=NA,alpha=0.8) +
    geom_dotplot(binaxis="y",stackdir="center",alpha=0.7,position=position_dodge(0.8)) +
    labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("\u03B4"^{bold("13")},"CO"[bold("2")]," (\u2030)"))),title = expression(bold("Leaf Substrate")),fill=expression(bold("Tree Species"))) +
    ylim(-20,100) +
    scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
    scale_fill_viridis(discrete=TRUE,labels=c("G. meintha","P. macroloba")) +
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
    geom_boxplot(lwd=1.5,outlier.shape=NA,alpha=0.8) +
    geom_dotplot(binaxis="y",stackdir="center",alpha=0.7,position=position_dodge(0.8)) +
    labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("\u03B4"^{bold("13")},"CO"[bold("2")]," (\u2030)"))),title = expression(bold("Starch Substrate")),fill=expression(bold("Tree Species"))) +
    ylim(-20,260) +
    scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
    scale_fill_viridis(discrete=TRUE,labels=c("G. meintha","P. macroloba")) +
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

  ## Separate leaf and starch datasets
  leaf_PLFA<-filter(data,Substrate_Type=="L")
  starch_PLFA<-filter(data,Substrate_Type=="S")  

  ## Leaf PLFA figure
  leaf_d13_PLFA<-ggplot(leaf_PLFA,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
    geom_hline(yintercept=-25,lwd=1.5,linetype="dashed") +
    geom_boxplot(lwd=1) +
    ylab(label=expression(bold(paste("\u03B4"^{bold("13")}, "C PLFA"," (\u2030)")))) +
    scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
    scale_fill_viridis(discrete=TRUE,labels=c("AMF","Bacteria","Fungi","All Microbes")) +
    xlab(label="Rhizosphere Manipulation") +
    ggtitle("Leaf Substrate") +
    labs(fill="Microbial Group") +
    facet_wrap(~Tree,labeller=labeller(Tree=tree_labs)) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          axis.text.y=element_text(colour="black",size=10),
          axis.text.x=element_text(colour="black",size=14),
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
    geom_hline(yintercept=-25,lwd=1.5,linetype="dashed") +
    geom_boxplot(lwd=1) +
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
          axis.text.x=element_text(colour="black",size=14),
          axis.title=element_text(size=14,face="bold"),
          plot.title=element_text(size=14,face="bold"),
          panel.border=element_rect(fill=NA,colour="black",size=1.5),
          panel.background=element_rect(fill=NA),
          strip.text.x=element_text(size=14,face="bold.italic"),
          legend.position="none")
  
  ## Panneled figure
  grid.arrange(leaf_d13_PLFA,starch_d13_PLFA,ncol=2)

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
  
  ## Respiration figure
  resp_fig<-ggplot(data=data,aes(x=Exclusion,y=Mean_Respiration)) +
    geom_boxplot(lwd=1,outlier.shape=NA) +
    geom_dotplot(binaxis="y",stackdir="center",alpha=0.7) +
    labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("Respiration Rate (",mu,"g C",O[2]," kg soi",l^-1,h^-1,")"))),title = expression(bold("Soil Respiration"))) +
    annotate("text", x = 1, y = 7, label = "A", size=8, color="red") +
    annotate("text", x = 2, y = 7, label = "A", size=8, color="red") +
    annotate("text", x = 3, y = 7, label = "B", size=8, color="red") +
    scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          axis.text.y=element_text(colour="black",size=10),
          axis.text.x=element_text(colour="black",size=14),
          axis.title=element_text(size=14,face="bold"),
          panel.border=element_rect(fill=NA,colour="black",size=1.5),
          panel.background=element_rect(fill=NA))
  
  ## ANOVA for soil respiration
  resp_mod<-lm(data=data,Mean_Respiration~Exclusion)
  Anova(resp_mod,type="II") # Exclusion highly significant, p = 0.006)
  TukeyHSD(aov(resp_mod)) # 3 is different from 1 and 2
  

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
  
  ### Take mean of narrow groups to yield FUNGI, BACTERIA, ALL_MICROBE groups
  for (i in 1:nrow(data_wide)) {
    data_wide$BACTERIA[i]<-sum(c(data_wide$bacteria[i],data_wide$`gram+bacteria`[i],data_wide$`gram-bacteria`[i]))
    data_wide$FUNGI[i]<-sum(c(data_wide$AMF[i],data_wide$fungi[i]))
    data_wide$ALL_MICROBES[i]<-sum(c(data_wide$AMF[i],data_wide$bacteria[i],data_wide$fungi[i],data_wide$general[i],data_wide$`gram-bacteria`[i],data_wide$`gram+bacteria`[i],data_wide$protozoa[i]))
  }
  
  data_wide<-data_wide[,c(1,2,3,4,11,12,13)] # Trim down to broader groups
  
  data<-melt(data_wide,value.name="umol_FAME",id.vars=c("Sample_Name","Tree","Rhizosphere_Manipulation"))
  names(data)[4]<-"Microbial_Group"
  
  ## ANOVA for general PLFA by rhizosphere manipulation and tree type
  PLFA_mod<-lm(data=data,umol_FAME~Rhizosphere_Manipulation*Microbial_Group+Tree)
  Anova(PLFA_mod,type="II") # Rhizospehre manipulation not significant, microbial group significant, tree significant, rhizosphere*microbial group not significant
  
  TukeyHSD(aov(PLFA_mod))
  HSD.test(y=aov(PLFA_mod),group=TRUE)
  
  