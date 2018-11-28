## Costa Rica 2016 13C PLFA
## Kenna Rewcastle
## 02/09/18

library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)

setwd("~/Desktop")
dat1<-read.csv("13C_PLFA_master.csv")
dat<-data.frame(dat1) # At some point, will need to filter out rows that don't have a compound name.


#### No need to re-run this part.
samp<-dat$Identifier.1
d13<-dat$d.13C.12C

d13_dat<-data.frame(samp,d13)
mean_dat<-ddply(d13_dat, "samp", numcolwise(mean))

check<-d13_dat[1:95,2] # mean_dat correctly averages d13 values for each sample across all lipids!
mean(check)

write.csv(mean_dat,file="mean_dat.csv")
####


d13_dat1<-read.csv("mean_dat.csv")
d13<-d13_dat1$d13
label<-d13_dat1$label_type
tree<-d13_dat1$tree.type
core<-as.factor(d13_dat1$core_type)

mod1<-lm(d13~core)
anova(mod1) # core has significant effect!!

d13_core<-ggplot(data=d13_dat1,aes(x=core,y=d13,fill=tree.type)) +
  geom_boxplot(lwd=1.5) +
  geom_hline(yintercept=-25, linetype="dashed",lwd=1.5) +
  labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("\u03B4"^{bold("13")}, "C (\u2030)"))),fill="Tree Species") +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  scale_fill_discrete(labels=c("G"="Goeth","P"="Penta")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

mod2<-lm(d13~label)
anova(mod2) # significant

d13_label<-ggplot(data=d13_dat1,aes(x=label,y=d13)) +
  geom_boxplot(lwd=1.5) +
  geom_hline(yintercept=-25, linetype="dashed",lwd=1.5) +
  labs(x=expression(bold("Substrate Added")),y=expression(bold(paste("\u03B4"^{bold("13")}, "C (\u2030)")))) +
  scale_x_discrete(labels=c("L"="Leaf Material","S"="Starch")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

mod3<-lm(d13~tree)
anova(mod3) # not significant

d13_tree<-ggplot(data=d13_dat1,aes(x=tree,y=d13)) +
  geom_boxplot(lwd=1.5) +
  geom_hline(yintercept=-25, linetype="dashed",lwd=1.5) +
  labs(x=expression(bold("Tree Type")),y=expression(bold(paste("\u03B4"^{bold("13")}, "C (\u2030)")))) +
  scale_x_discrete(labels=c("P"="Pentaclethra (Legume)","G"="Goethalsia")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

###### Separate everything out by substrate type
d13leaf<-filter(d13_dat1,label=="L")
d13starch<-filter(d13_dat1,label=="S")

d13L<-d13leaf$d13
d13S<-d13starch$d13
coreL<-as.factor(d13leaf$core_type)
coreS<-as.factor(d13starch$core_type)
treeL<-d13leaf$tree.type
treeS<-d13starch$tree.type

# effect of rhizosphere manipulation on leaf decomp

mod4<-lm(d13L~coreL*treeL) # super super significant
anovaMod4<-aov(mod4)
anova(mod4)
TukeyHSD(anovaMod4) # mesh 3 is signficantly different from the other two, which are identical (root presence is the important driver here???)

leaf_d13_core<-ggplot(data=d13leaf,aes(x=coreL,y=d13L,fill=tree.type)) +
  geom_boxplot(lwd=1.5) +
  geom_hline(yintercept=-25, linetype="dashed",lwd=1.5) +
  labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("\u03B4"^{bold("13")}, "C  in Microbial PLFA (\u2030)"))),title = expression(bold("Leaf Substrate")),fill="Tree Species") +
  #annotate("text", x = 1, y = 27, label = "A", size=8, color="red") +
  #annotate("text", x = 2, y = 27, label = "A", size=8, color="red") +
  #annotate("text", x = 3, y = 27, label = "B", size=8, color="red") +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  scale_fill_discrete(labels=c("G"="Goeth","P"="Penta")) +
  ylim(-30,27) +
  guides(fill=FALSE) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

# effect of rhizosphere manipulation on starch decomp

mod5<-lm(d13S~coreS*treeS) # Not significant!!!
anova(mod5)
anovaMod5<-aov(mod5)
TukeyHSD(anovaMod5) # no differences among mesh

starch_d13_core<-ggplot(data=d13starch,aes(x=coreS,y=d13S,fill=tree.type)) +
  geom_boxplot(lwd=1.5) +
  geom_hline(yintercept=-25, linetype="dashed",lwd=1.5) +
  labs(x=expression(bold("Rhizosphere Manipulation")),y=expression(bold(paste("\u03B4"^{bold("13")}, "C in Microbial PLFA (\u2030)"))), title = expression(bold("Starch Substrate")),fill="Tree Species") +
  scale_x_discrete(labels=c("1"="-R-M","2"="-R+M","3"="+R+M")) +
  scale_fill_discrete(labels=c("G"="Goeth","P"="Penta")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))


grid.arrange(leaf_d13_core, starch_d13_core, nrow = 1)

# effect of tree type on leaf decomp

mod6<-lm(d13L~treeL)
anova(mod6) # not significant

leaf_d13_tree<-ggplot(data=d13leaf,aes(x=treeL,y=d13L)) +
  geom_boxplot(lwd=1.5) +
  geom_hline(yintercept=-25, linetype="dashed",lwd=1.5) +
  labs(x=expression(bold("Tree Type")),y=expression(bold(paste("\u03B4"^{bold("13")}, "C (\u2030)")))) +
  scale_x_discrete(labels=c("P"="Pentaclethra (Legume)","G"="Goethalsia")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))

# effect of tree type on starch decomp
mod7<-lm(d13S~treeS)
anova(mod7) # not significant

starch_d13_tree<-ggplot(data=d13starch,aes(x=treeS,y=d13S)) +
  geom_boxplot(lwd=1.5) +
  geom_hline(yintercept=-25, linetype="dashed",lwd=1.5) +
  labs(x=expression(bold("Tree Type")),y=expression(bold(paste("\u03B4"^{bold("13")}, "C (\u2030)")))) +
  scale_x_discrete(labels=c("P"="Pentaclethra (Legume)","G"="Goethalsia")) +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.text.y=element_text(colour="black",size=10),
        axis.text.x=element_text(colour="black",size=14),
        axis.title=element_text(size=14,face="bold"),
        panel.border=element_rect(fill=NA,colour="black",size=1.5),
        panel.background=element_rect(fill=NA))
  