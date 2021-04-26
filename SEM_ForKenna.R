arc<- read.csv("ARC_Diversity_Function.csv", header=TRUE)
names(arc)
library(lavaan)
library(qgraph)
library(semPlot)


#One way scale variables for SEM

zDis<-scale(arc$Fdis)
zDiv<-scale(arc$Fdiv)
zEve<-scale(arc$Feve)
zRic<-scale(arc$Fric)
zMPD<-scale(arc$MPD)
zFaith<-scale(arc$Faith.s)
zSR<-scale(arc$SR)
zH<-scale(arc$H)
zD<- scale(arc$D)
zNEE<-scale(arc$NEE)
zTB<-scale(arc$TotalBiomass)
zVWC<- scale(arc$VWC)

df= data.frame(zDis, zDiv, zEve, zRic, zMPD, zFaith, zSR, zH, zD, zNEE, zTB, zVWC)

model<- '
#measurement model (latent variables created here)

FD =~ 1 * zEve + zDis
EF =~ 1 * zNEE + zTB

#regressions
EF ~ FD + zH + zFaith + zRic


#residual correlations if you have any

'
#fit SEM

fit<- sem(model, data = df)
summary(fit, standardized = TRUE, rsq=T)

#plot results using semPaths in qgraph
semPaths(fit, "std", edge.label.cex = 0.5, curvePivot = TRUE, layout = "tree2", posCol = "navy", negCol = "peru", colFactor = .1)

#Example using your own data/variable names

model<- '
#measurement model (latent variables created here)

MicroAct =~ 1 * CO2 + Enzymes

#regressions
MicroAct ~ MicroNMDS + MicroCommStress + rootbiomass + mycoBiomass
MicroAct ~ mycoBiomass
MicroAct ~ rootbiomass
MicroCommStress ~ rootbiomass
MicroCommStress ~ mycoBiomass
MicroNMDS ~ rootbiomass

#residual correlations if you have any
mycoBiomass ~ rootbiomass
'


##Another way to scale - Example 
dat<- read.csv("Dimensions_RMBL_Core.csv", header=TRUE)
names(dat)
subset<- select(dat, GPP, NEE)

select(NEE, MPD, PD, SR, FDis, summerprecip, summertemp, AET, Elevation, Total_Bio, LAI) 

df2= dat %>% 
select(NEE, MPD, PD, SR, FDis, summerprecip, summertemp, AET, Elevation, Total_Bio, LAI) %>%
mutate(logSR=log(SR), logMPD=log(MPD), logprecip=log(summerprecip), logtemp=log(summertemp), logElevation=log(Elevation), logAET=log(AET), logPD=log(PD), logFDis=log(FDis), logBio=log(Total_Bio), NEEt=log(-1*NEE), logLAI=log(LAI)) %>%
scale() %>%
data.frame()
df2 <- data.frame(df2, dat$Site, dat$Quadrat, dat$Elevation)
#df2$fElevation <- factor(df2$Elevation) 
head(df2)


