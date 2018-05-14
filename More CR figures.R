endoDat<-read.csv("CR_2016_fungal_colonization_DOE_cores.csv")
names(endoDat)
endoDat$Core<-as.factor(endoDat$Core)
endoDat<-filter(endoDat,Core!="NA")

ggplot(endoDat,aes(x=Core,y=per_AMF)) +
  geom_boxplot()

mod<-lm(endoDat$per_AMF~endoDat$Core)
anova(mod)
anovaMod<-aov(mod)
TukeyHSD(anovaMod)

names(endoDat)

DatDat<-read.csv("Costa Rica Master sheet annotated.csv")
names(DatDat)
hyphal_scaled<-DatDat$Dry_root_wt*DatDat$per_AMF
newDat<-data.frame(DatDat,hyphal_scaled)
newDat$Exclusion<-as.factor(newDat$Exclusion)
newDat<-filter(newDat,Exclusion!="NA")

ggplot(newDat,aes(x=Exclusion,y=hyphal_scaled)) + # y-axis
  geom_boxplot()
