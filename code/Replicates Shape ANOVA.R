library(geomorph)
setwd("~/MCGILL/Bluegill Sunfish Data/2018 Photos/BLUEGILL/2018 BLG BODY PHOTOS/Replicates")
myData_tpsrep <- readland.tps("Total Replicates Coords.TPS ",specID="imageID")
myData_tpsrep

### 1. Make Data Frame
gdf.fish<-geomorph.data.frame(shape=myData_tpsrep)

myData.superimp=gpagen(myData_tpsrep, ProcD = TRUE, Proj = TRUE)  
plot(myData.superimp)
corePCA=plotTangentSpace(myData.superimp$coords) ###Page57
View(corePCA)
summary(corePCA)


### Classifiers
classifiers=read.table("Replicates Identifiers.txt",sep = "\t", header=TRUE)
replic<-classifiers$replicateNo
plotTangentSpace(myData.superimp$coords, groups=as.factor(paste(replic,classifiers$lakeID)), legend = TRUE)

gdf.fish$replic<-as.factor(classifiers$replicateNo)
gdf.fish=geomorph.data.frame(shape=myData_tpsrep,replic=classifiers$replicateNo,fishID=classifiers$imageID,lake=classifiers$lakeID,DOC=classifiers$lakeDOC,cSize=myData.superimp$Csize)

#ANOVA
fishrep.anova<-procD.lm(shape ~ fishID/replic, data=gdf.fish) 
summary(fishrep.anova)
fishrep.anova$coefficients

###             Df     SS      MS     Rsq        F       Z Pr(>SS)   
#fishID         59 602102 10205.1 0.97561 187.1758 10.3693   0.001 **
# fishID:replic  60   8509   141.8 0.01379   2.6012  5.3186   0.001 **
#Residuals     120   6543    54.5 0.01060                            
#Total         239 617154                    


# Percent repeatability. Grant and Madlen did this. 
part1<-(fishrep.anova$aov.table[1,3]-fishrep.anova$aov.table[2,3])/4

part1/(fishrep.anova$aov.table[2,3]+part1)





