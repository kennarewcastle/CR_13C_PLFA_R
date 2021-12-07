#### Current Costa Rica Figures
#### KER
#### Created on: 07 December 2021


# Load packages -----------------------------------------------------------

library(dplyr)
library(ggplot2)
library(gridExtra)
library(reshape2)

# delta13C-CO2 ------------------------------------------------------------

  ## Read in data
  dat<-read.csv("Costa Rica Master sheet annotated.csv")
  dat<-filter(dat,dat$Exclusion!="NA")
  dat$Exclusion<-as.factor(dat$Exclusion)

  ## Separate starch and leaf data
  starch_CO2<-filter(dat,dat$Isotope_label=="S") 
  leaf_CO2<-filter(dat,dat$Isotope_label=="L")

  ## Leaf d13 figure
  leaf_d13<-ggplot(data=leaf_CO2,aes(x=Exclusion,y=d13_d9)) +
    geom_boxplot(lwd=1) +
    geom_point(alpha=0.7,size=5) +
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

  ## 13C CO2 leaf anova and TukeyHSD
  leafMod<-lm(leaf_CO2$d13_d9~leaf_CO2$Exclusion) # super super significant p < 0.001
  anova(leafMod) # Significant rhizosphere manipulation
  TukeyHSD(aov(leafMod)) # -R-M is different from the other two exclusion types

  ## Starch d13 figure
  starch_d13<-ggplot(data=starch_CO2,aes(x=Exclusion,y=d13_d2)) +
    geom_boxplot(lwd=1) +
    geom_point(size=5, alpha=0.7) +
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

  ## 13C CO2 starch anova and TukeyHSD
  starchMod<-lm(starch_CO2$d13_d2~starch_CO2$Exclusion)
  anova(starchMod) # Only significant differences are exclusion
  TukeyHSD(aov(starchMod))

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
  data_trim<-data.frame("Sample_Name"=data[,1],"Rhizosphere_Manipulation"=data[,3],"Substrate_Type"=data[,4],"Microbial_Group"=data[,14],"d13"=data[,9])
  
  data_wide<-dcast(data_trim,Sample_Name+Rhizosphere_Manipulation+Substrate_Type ~ Microbial_Group, value.var="d13", fun.aggregate=mean, na.rm=TRUE)
  data_wide[is.nan(data_wide)] <- NA
  
  data_wide$BACTERIA<-mean(c(data_wide$bacteria,data_wide$`gram+bacteria`,data_wide$`gram-bacteria`),na.rm=TRUE)
  data_wide$FUNGI<-mean(c(data_wide$AMF,data_wide$fungi),na.rm=TRUE)
  data_wide$ALL_MICROBES<-mean(c(data_wide$AMF,data_wide$bacteria,data_wide$fungi,data_wide$general,data_wide$`gram-bacteria`,data_wide$`gram+bacteria`),na.rm=TRUE)
  data_wide<-data_wide[,c(1,2,3,4,11,12,13)]
  
  data<-melt(data_wide,value.name="d13",id.vars=c("Sample_Name","Rhizosphere_Manipulation","Substrate_Type"))
  names(data)[4]<-"Microbial_Group"

  ## Separate leaf and starch datasets
  leaf_PLFA<-filter(data,Substrate_Type=="L")
  starch_PLFA<-filter(data,Substrate_Type=="S")  

  ## Leaf PLFA figure
  leaf_d13_PLFA<-ggplot(leaf_PLFA,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
    geom_boxplot(lwd=1) +
    scale_fill_brewer(palette = "Set1") +
    ylab(label=expression(bold(paste("\u03B4"^{bold("13")}, "C PLFA"," (\u2030)")))) +
    scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
    xlab(label="Rhizosphere Manipulation") +
    ggtitle("Leaf Substrate") +
    labs(fill="Microbial Group") +
    #ylim(-25,75) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          axis.text.y=element_text(colour="black",size=10),
          axis.text.x=element_text(colour="black",size=14),
          axis.title=element_text(size=14,face="bold"),
          plot.title=element_text(size=14,face="bold"),
          panel.border=element_rect(fill=NA,colour="black",size=1.5),
          panel.background=element_rect(fill=NA),
          legend.position="none")
  
  ## Starch PLFA figure
  starch_d13_PLFA<-ggplot(starch_PLFA,aes(x=Rhizosphere_Manipulation,y=d13,fill=Microbial_Group)) +
    geom_boxplot(lwd=1) +
    scale_fill_brewer(palette = "Set1") +
    ylab(label=expression(bold(paste("\u03B4"^{bold("13")}, "C PLFA"," (\u2030)")))) +
    scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
    xlab(label="Rhizosphere Manipulation") +
    ggtitle("Starch Substrate") +
    labs(fill="Microbial Group") +
    #ylim(-50,275) +
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
          axis.text.y=element_text(colour="black",size=10),
          axis.text.x=element_text(colour="black",size=14),
          axis.title=element_text(size=14,face="bold"),
          plot.title=element_text(size=14,face="bold"),
          panel.border=element_rect(fill=NA,colour="black",size=1.5),
          panel.background=element_rect(fill=NA),
          legend.position="none")
  
  ## Panneled figure
  grid.arrange(leaf_d13_PLFA,starch_d13_PLFA,ncol=2)
  