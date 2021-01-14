# This script creates SVG files showing morphometric landmarks.

# Load packages -----------------------------------------------------------
library(geomorph)
library(shapes)
library(Morpho)
library(StereoMorph)
library(phangorn)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(here)

# Load data ---------------------------------------------------------------
myData_tps <- readland.tps(here("data", "FULL_2018_TPS_FILE_UPDATED_09-25-19.TPS"), specID = "imageID")
mylinks <- read.table(here("data", "Full_body_links.txt"))
identifiers <- read.table(here("data", "Identifiers_Update_2020.txt"), sep = "\t", header=TRUE)

dimnames(myData_tps)[[3]] <- identifiers$imageID

dimnames(myData_tps)[[1]] <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19") 

GPA.fish <- gpagen(myData_tps, ProcD = TRUE, Proj = TRUE) 

gdf.fish <- geomorph.data.frame(shape=GPA.fish$coords,
                                DOCrange=identifiers$DOCrange,
                                Lake=identifiers$lakeID,
                                Sex=identifiers$Sex, 
                                DOC=identifiers$lakeDOC, 
                                DOCcat=identifiers$DOC,
                                captureMethod=identifiers$captureMethod, 
                                cSize=GPA.fish$Csize, 
                                basin=identifiers$Basin)

corePCA <- gm.prcomp(A = gdf.fish$shape, # XXX come back to this!
                            groups=as.factor(gdf.fish$DOC),
                            legend = TRUE)

 #XXX START HERE
PCscores=corePCA$pc.scores[,1:2]
idx=which(gdf.fish$Lake=="Bay"); Bay=PCscores[idx,]
idx=which(gdf.fish$Lake=="Birch"); Birch=PCscores[idx,]
idx=which(gdf.fish$Lake=="Crampton"); Crampton=PCscores[idx,]
idx=which(gdf.fish$Lake=="Found"); Found=PCscores[idx,]
idx=which(gdf.fish$Lake=="Hummingbird"); Hummingbird=PCscores[idx,]
idx=which(gdf.fish$Lake=="Little_Crooked"); Little_Crooked=PCscores[idx,]
idx=which(gdf.fish$Lake=="Lost"); Lost=PCscores[idx,]
idx=which(gdf.fish$Lake=="McCullough"); McCullough=PCscores[idx,]
idx=which(gdf.fish$Lake=="Muskellunge"); Muskellunge=PCscores[idx,]
idx=which(gdf.fish$Lake=="Oxbow"); Oxbow=PCscores[idx,]
idx=which(gdf.fish$Lake=="Papoose"); Papoose=PCscores[idx,]
idx=which(gdf.fish$Lake=="Red_Bass"); Red_Bass=PCscores[idx,]
idx=which(gdf.fish$Lake=="Squaw"); Squaw=PCscores[idx,]
idx=which(gdf.fish$Lake=="Towanda"); Towanda=PCscores[idx,]
plotRefToTarget(corePCA$pc.shapes$PC1min, corePCA$pc.shapes$PC1max,method="points", links = mylinks, gridPars = gridPar(tar.pt.bg = "black",tar.link.col = "black",tar.link.lwd = 3, tar.pt.size = 1, pt.size = 1, pt.bg = "gray", link.lwd = 3), mag=1, useRefPts = TRUE)
plotRefToTarget(corePCA$pc.shapes$PC2min, corePCA$pc.shapes$PC2max,method="points", links = mylinks, gridPars = gridPar(tar.pt.bg = "black",tar.link.col = "black",tar.link.lwd = 3, tar.pt.size = 1, pt.size = 1, pt.bg = "gray", link.lwd = 3), mag=1, useRefPts = TRUE)

library(StereoMorph)
lm_array<-myData_tps
gpa_array <- gpagen(myData_tps, ProcD = TRUE, Proj = TRUE)$coords ### takes just the coords

# Convert array to matrix for PCA
gpa_mat <- t(apply(gpa_array, 3, function(y) matrix(t(y), 1)))

# Perform non-phylogenetic PCA
resEig <- eigen(cov(gpa_mat))

# Get PC scores, changed sign of rotation matrix (resEig$vectors)

scores <- gpa_mat %*% -(resEig$vectors)
####shapes switched, need to reverse signs to fix

# Get percent variance explained along each axis
per_var <- (resEig$values / sum(resEig$values))*100

plot_fish_lateral <- function(xy, coor, size=1, col="black"){
  
  # If 3D, rotate points about x-axis using 3D rotation matrix
  if(ncol(coor) == 3){
    coor <- coor %*% matrix(c(1,0,0, 0,cos(-pi/2),sin(-pi/2), 
                              0,-sin(-pi/2),cos(-pi/2)), nrow=3, ncol=3)
  }
  
  # Get just x,y coordinates (orthographic projection into xy-plane)
  coor <- coor[, 1:2]
  
  # Get plot aspect ratio
  w <- par("pin")[1]/diff(par("usr")[1:2])
  h <- par("pin")[2]/diff(par("usr")[3:4])
  asp <- w/h
  
  # Correct for plot aspect ratio not necessarily being 1:1
  coor[, 1] <- coor[, 1] * (1/asp)
  
  # Scale points and place back in position
  coor <- coor*size
  
  # Center about zero based on range of coordinates
  coor <- coor - matrix(colMeans(apply(coor, 2, range)), 
                        nrow=nrow(coor), ncol=ncol(coor), byrow=TRUE)
  
  # Move shape to PC score
  coor <- coor + matrix(xy, nrow(coor), ncol(coor), byrow=TRUE)
  
  # Set order in which to draw points to create polygon
  polygon_order=c(1,3:12,16,1) # name landmarks in the order you want them connected
  # Create filled polygon
  polygon(coor[polygon_order, ], col=col, border=col)
}
# Set PCs to plot
pcs <- 1:2

# prepare input data
gdf.fish=geomorph.data.frame(shape=GPA.fish$coords,DOCrange=identifiers$DOCrange,
                             Lake=identifiers$lakeID,Sex=identifiers$Sex, DOC=identifiers$lakeDOC, DOCcat=identifiers$DOC,
                             captureMethod=identifiers$captureMethod, cSize=GPA.fish$Csize, basin=identifiers$Basin)

PCscores=corePCA$pc.scores[,1:2] # which groups do you want to plot? put basin, when basin
#idx=which(gdf.fish$DOCcat=="low"); Low_DOC =PCscores[idx,]
#idx=which(gdf.fish$DOCcat=="high"); High_DOC=PCscores[idx,]

#
idx=which(gdf.fish$Lake=="Bay"); Bay=PCscores[idx,]
idx=which(gdf.fish$Lake=="Birch"); Birch=PCscores[idx,]
idx=which(gdf.fish$Lake=="Crampton"); Crampton=PCscores[idx,]
idx=which(gdf.fish$Lake=="Found"); Found=PCscores[idx,]
idx=which(gdf.fish$Lake=="Hummingbird"); Hummingbird=PCscores[idx,]
idx=which(gdf.fish$Lake=="Little_Crooked"); Little_Crooked=PCscores[idx,]
idx=which(gdf.fish$Lake=="Lost"); Lost=PCscores[idx,]
idx=which(gdf.fish$Lake=="McCullough"); McCullough=PCscores[idx,]
idx=which(gdf.fish$Lake=="Muskellunge"); Muskellunge=PCscores[idx,]
idx=which(gdf.fish$Lake=="Oxbow"); Oxbow=PCscores[idx,]
idx=which(gdf.fish$Lake=="Papoose"); Papoose=PCscores[idx,]
idx=which(gdf.fish$Lake=="Red_Bass"); Red_Bass=PCscores[idx,]
idx=which(gdf.fish$Lake=="Squaw"); Squaw=PCscores[idx,]
idx=which(gdf.fish$Lake=="Towanda"); Towanda=PCscores[idx,]

# Create plot box with axes and axis labels
plot(scores[, pcs], type="n", main="Backtransform morphospace",
     xlab=paste0("PC", pcs[1], " (", round(per_var[pcs[1]]), "%)"),
     ylab=paste0("PC", pcs[2], " (", round(per_var[pcs[2]]), "%)"))

# Plot backtransform shapes, changed sign of rotation matrix (resEig$vectors) 
btShapes(scores=scores, vectors=-(resEig$vectors), fcn=plot_fish_lateral, 
         pcs=pcs, n=c(4,4), m=dim(lm_array)[2], row.names=dimnames(lm_array)[[1]], 
         pc.margin=c(0.06,0.05), size=0.038, col=gray(0.7))
#points(Low_DOC [,1], Low_DOC[,2], col="#41C6EC", pch=19,cex=2)
#points(High_DOC [,1] , High_DOC[,2], col="#A0674E", pch=19,cex=2)
#legend(0.05,0.04, legend=c("Low_DOC","High_DOC"), pch=19, col=c("#41C6EC","#A0674E"))
points(Bay[,1], Bay[,2], col="#6666FF", pch=19,cex=2)
points(Birch[,1], Birch[,2], col="#3333FF", pch=19,cex=2)
points(Crampton[,1], Crampton[,2], col="#33FFFF", pch=19,cex=2)
points(Found[,1], Found[,2], col="#33CCFF", pch=19,cex=2)
points(Hummingbird[,1], Hummingbird[,2], col="#000033", pch=19,cex=2)
points(Little_Crooked[,1], Little_Crooked[,2], col="#99FFFF", pch=19,cex=2)
points(Lost[,1], Lost[,2], col="#CCFFFF", pch=19,cex=2)
points(McCullough[,1], McCullough[,2], col="#333399", pch=19,cex=2)
points(Muskellunge[,1], Muskellunge[,2], col="#6699FF", pch=19,cex=2)
points(Oxbow[,1], Oxbow[,2], col="#3333CC", pch=19,cex=2)
points(Papoose[,1], Papoose[,2], col="#0099FF", pch=19,cex=2)
points(Red_Bass[,1], Red_Bass[,2], col="#003366", pch=19,cex=2)
points(Squaw[,1], Squaw[,2], col="#333366", pch=19,cex=2)
points(Towanda[,1], Towanda[,2], col="#99CCFF", pch=19,cex=2)
legend("topright", legend=c("Bay","Birch","Crampton", "Found", "Hummingbird", "Little_Crooked", "Lost", "McCullough", 
                            "Muskellunge", "Oxbow", "Papoose", "Red_Bass", "Squaw", "Towanda" ), 
       pch=19, col=c("#6666FF","#3333FF","#33FFFF","#33CCFF","#000033","#99FFFF", "#CCFFFF", "#333399", 
                     "#6699FF", "#3333CC", "#0099FF", "#003366", "#333366", "#99CCFF", cex=0.30 ))

# Add convex hull polygons to the PCA plot, PC1 vs PC2:
colour = c("#6666FF","#3333FF","#33FFFF","#33CCFF","#000033","#99FFFF", "#CCFFFF", "#333399", 
           "#6699FF", "#3333CC", "#0099FF", "#003366", "#333366", "#99CCFF")
for(j in 1:nlevels(as.factor(gdf.fish$DOC))) {
  # Get edge points (used to plot convex hull):
  edge_points <- rownames(PCscores[which(gdf.fish$DOC == levels(as.factor(gdf.fish$DOC))[j]),])[
    chull(PCscores[which(gdf.fish$DOC == levels(as.factor(gdf.fish$DOC))[j]), c(1,2)])]
  # Plot convex hull as polygon:
  polygon(PCscores[edge_points, c(1,2)], col = adjustcolor(colour[j],
                                                           alpha.f = 0.3) , border = colour[j])
} # alpha gives the degree of transparency of the polygon

###Mean Shapes Per Lake
lm_array<-myData_tps
gpa_array <- gpagen(myData_tps, ProcD = TRUE, Proj = TRUE)$coords ### takes just the coords

x <- two.d.array(gpa_array) # the following part is necessary to calculate mean shapes per lake
gpa_df=geomorph.data.frame(shape=gpa_array, Lake=identifiers$lakeID)
p <- dim(gpa_array)[1] # the number of landmarks
k <- dim(gpa_array)[2] # the dimensions of the coordinates
Y <- array(NA, dim=c(p,k,length(levels(gpa_df$Lake)))) #new empty array to fill
dimnames(Y)[[3]] <- levels(gpa_df$Lake) # set group levels as new names

for (i in 1: length(levels(gpa_df$Lake))){
  grp <- x[which(gpa_df$Lake==levels(gpa_df$Lake)[i]),]
  foo <- arrayspecs(grp ,p,k)
  Y[,,i] <- mshape(foo) # place into the new 3D array
}
gpa_array <- Y

# Convert array to matrix for PCA
gpa_mat <- t(apply(gpa_array, 3, function(y) matrix(t(y), 1)))

# Perform non-phylogenetic PCA
resEig <- eigen(cov(gpa_mat))

# Get PC scores, changed sign of rotation matrix (resEig$vectors)

scores <- gpa_mat %*% -(resEig$vectors)
####shapes switched need to reverse signs to fix

# Get percent variance explained along each axis
per_var <- (resEig$values / sum(resEig$values))*100

plot_fish_lateral <- function(xy, coor, size=1, col="black"){
  
  # If 3D, rotate points about x-axis using 3D rotation matrix
  if(ncol(coor) == 3){
    coor <- coor %*% matrix(c(1,0,0, 0,cos(-pi/2),sin(-pi/2), 
                              0,-sin(-pi/2),cos(-pi/2)), nrow=3, ncol=3)
  }
  
  # Get just x,y coordinates (orthographic projection into xy-plane)
  coor <- coor[, 1:2]
  
  # Get plot aspect ratio
  w <- par("pin")[1]/diff(par("usr")[1:2])
  h <- par("pin")[2]/diff(par("usr")[3:4])
  asp <- w/h
  
  # Correct for plot aspect ratio not necessarily being 1:1
  coor[, 1] <- coor[, 1] * (1/asp)
  
  # Scale points and place back in position
  coor <- coor*size
  
  # Center about zero based on range of coordinates
  coor <- coor - matrix(colMeans(apply(coor, 2, range)), 
                        nrow=nrow(coor), ncol=ncol(coor), byrow=TRUE)
  
  # Move shape to PC score
  coor <- coor + matrix(xy, nrow(coor), ncol(coor), byrow=TRUE)
  
  # Set order in which to draw points to create polygon
  polygon_order=c(1,3:12,16,1) # name landmarks in the order you want them connected
  # Create filled polygon
  polygon(coor[polygon_order, ], col=col, border=col)
}
# Set PCs to plot
pcs <- 1:2

# Create plot box with axes and axis labels
plot(scores[, pcs], type="n", main="Backtransform morphospace",
     xlab=paste0("PC", pcs[1], " (", round(per_var[pcs[1]]), "%)"),
     ylab=paste0("PC", pcs[2], " (", round(per_var[pcs[2]]), "%)"))


# Plot backtransform shapes, changed sign of rotation matrix (resEig$vectors) 
btShapes(scores=scores, vectors=-(resEig$vectors), fcn=plot_fish_lateral, 
         pcs=pcs, n=c(4,4), m=dim(lm_array)[2], row.names=dimnames(lm_array)[[1]], 
         pc.margin=c(0.06,0.06), size=0.018, col=gray(0.7))
points(scores [,1], scores[,2], col="#41C6EC", pch=19,cex=2)
text(scores [,1], scores [,2], labels=levels(gpa_df$Lake), cex= 0.7)

##Ref to Target Figures
plotRefToTarget(corePCA$pc.shapes$PC1min, corePCA$pc.shapes$PC1max,method="points", links = mylinks, gridPars = gridPar(tar.pt.bg = "black",tar.link.col = "black",tar.link.lwd = 3, tar.pt.size = 1, pt.size = 1, pt.bg = "gray", link.lwd = 3), mag=1, useRefPts = TRUE)
plotRefToTarget(corePCA$pc.shapes$PC2min, corePCA$pc.shapes$PC2max,method="points", links = mylinks, gridPars = gridPar(tar.pt.bg = "black",tar.link.col = "black",tar.link.lwd = 3, tar.pt.size = 1, pt.size = 1, pt.bg = "gray", link.lwd = 3), mag=1, useRefPts = TRUE)
----------------------------------------------

##RE-RUN STATS FOR SHAPE TO CHECK THEM
final=anova(procD.lm(shape ~ log(cSize) + DOC + Lake + 
                       log(cSize)*DOC ,SS.type = "II",data = gdf.fish),
            error= c("Residuals","Residuals","Lake","Residuals"))
summary(final)
finalWBasin=anova(procD.lm(shape ~ log(cSize) + DOC + basin/Lake +
                       log(cSize)*DOC ,SS.type = "II",data = gdf.fish),
            error= c("Residuals","Residuals","basin","basin:Lake","Residuals"))
summary(finalWBasin)
finalWBasin2=anova(procD.lm(shape ~ log(cSize) + DOC + basin/Lake +
                             log(cSize)*DOC ,SS.type = "II",data = gdf.fish),
                  error= c("Residuals","Residuals","Residuals","basin:Lake","Residuals"))
summary(finalWBasin2)
-----------------------------------------------------

###Univariate Data Check and Run

#Eyes
setwd("~/")
eyeWidthDat<-read.csv("eyewidthsFINAL.csv")

dfeye<-data.frame(eyewidth<-eyeWidthDat$eyeWidth,
                  fishIDeye<-eyeWidthDat$fishID,
                  lakeID<-eyeWidthDat$lakeID,
                  fishLength<-eyeWidthDat$fishStdLength,
                  DOClevel<-eyeWidthDat$DOClevel,
                  DOC<-eyeWidthDat$DOC,
                  eyewidth.ss<-eyeWidthDat$sizeStandardize,
                  lakeSS<-eyeWidthDat$lakeSS,basin<-eyeWidthDat$basin, 
                  area<-eyeWidthDat$area,maxDepth<-eyeWidthDat$max_Depth, fittedValues<-eyeWidthDat$fitted, fitted2<-eyeWidthDat$fitted2)
#Original Model 
EyeModel=lmer(log(eyewidth.ss) ~ 1 + log(DOC) + (1|lakeID), data=dfeye)
summary(EyeModel)
#Model with basin. Lake is nested within basin and set as random
EyeModelWBasin=lmer(log(eyewidth.ss) ~ 1 + log(DOC) + (1|basin/lakeID), data=dfeye)
summary(EyeModelWBasin)
#Try Model with intercepts varying by basin:lake and slopes varying by DOC in basin
EyeBasinNew=lmer(log(eyewidth.ss) ~ 1 + log(DOC) + basin + (1|basin:lakeID) + (log(DOC)|basin), data=dfeye)
summary(EyeBasinNew)
EyeBasinNew2=(lmer(log(eyewidth.ss) ~ (1|basin:lakeID) + (log(DOC)|basin), data=dfeye))
(EyeBasinNew2)



#get fitted values and add to data frame. Note only one value for each lake, so only 14 different fitted values. fitted(model)
fitted(EyeModel)
fitted(EyeModelWBasin)
##Same fitted values


ggplot(dfeye,aes(x = DOC, y =fittedValues, label= lakeID)) + 
  geom_point(aes(colour = lakeID)) + 
  geom_smooth(method = 'lm')
---------------------------------------------------------------------

#Gill Rakers
rakerData<- read.csv("Gill_Rakers_2018_Final.csv")
dfraker<-data.frame(captureMethod<-rakerData$capture_Method,lakeLss=rakerData$avgrakerlengthSS_bylake,
                    lakeSss=rakerData$avgrakerspaceSS_bylake,lakeCss=rakerData$rakercountSS_bylake,
                    rakercount=rakerData$total_RakerNum,DOClevel=rakerData$DOCbin,DOC=rakerData$lakeDOC,
                    Lake=rakerData$lakeID,fishLength=rakerData$fishSL,avgRakerLength=rakerData$avgRakerLength,
                    avgRakerSpace=rakerData$avgRakerSpace,count.ss=rakerData$SS.Count, length.ss=rakerData$SS.Length,
                    space.ss=rakerData$SS.Space,avgL2<-rakerData$avgL_4.7,avgS2<-rakerData$avgS_4.6,
                    avgL2_ss<-rakerData$avgL2_ss,avgS2_ss<-rakerData$avgS2_ss, basin.raker<-rakerData$basin, fitted_L<-rakerData$fitted_L,fitted_S<-rakerData$fitted_S,fitted_C<-rakerData$fitted_C)
#avgL2_ss is size standardizations for average length for rakers 4-7
rakerLModel=lmer(log10(avgL2_ss) ~ 1 + log10(DOC) + (1|Lake), data=dfraker)
summary(rakerLModel)
rakerSModel=lmer(log10(avgS2_ss) ~ 1 + log10(DOC) + (1|Lake), data=dfraker)
summary(rakerSModel)
#with Lake nested within basin, set as random
rakerLModelB=lmer(log10(avgL2_ss) ~ 1 + log10(DOC) + (1|basin.raker/Lake), data=dfraker)
summary(rakerLModelB)
rakerSModelB=lmer(log10(avgS2_ss) ~ 1 + log10(DOC) + (1|basin.raker/Lake), data=dfraker)
summary(rakerSModelB)
#NOTE: Higher p values when basin included


#use glm for count data
rakerCModel=glmer(rakercount ~ 1 + log10(DOC) + (1|Lake),family = poisson,nAGQ=0,control=glmerControl(optimizer = "nloptwrap"), data=dfraker)
anova(rakerCModel)
summary(rakerCModel)
rakerCModellog=glmer(log10(rakercount) ~ 1 + log10(DOC) + (1|Lake),family = poisson,nAGQ=0,control=glmerControl(optimizer = "nloptwrap"), data=dfraker)
anova(rakerCModellog)
summary(rakerCModellog)
#In RMarkdown, proper ANOVA table shows up when rakercount is logged, here it does not. I think this is because it is at the top of all of the "non-integers" but there is not enough room in the console to show it.
#Below are results from logged model in RMarkdown
#Fixed effects:
#              Estimate Std. Error z value Pr(>|z|)
#(Intercept)  0.06687    1.01802   0.066    0.948
#log10(DOC)  -0.01873    1.03944  -0.018    0.986
#No significant effect of DOC on count
#Results from rakerCModel (not logged count)
#            Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  2.46900    0.05549  44.497   <2e-16 ***
#  log10(DOC)  -0.04525    0.05696  -0.794    0.427 
rakerCModelLM=lmer(rakercount ~ 1 + log10(DOC) + (1|Lake), data=dfraker)
summary(rakerCModelLM)
#Just in case, ran lmer as well, still no significant effect of DOC on count

#with basin - lake nested within basin and set to random
rakerCModelB=glmer(rakercount ~ 1 + log10(DOC) + (1|basin.raker/Lake),family = poisson,nAGQ=0,control=glmerControl(optimizer = "nloptwrap"), data=dfraker)
anova(rakerCModelB)
summary(rakerCModelB)
#Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  2.46900    0.05549  44.497   <2e-16 ***
#  log10(DOC)  -0.04525    0.05696  -0.794    0.427 
#No Effect of DOC on raker number

#get fitted values for figures
fitted(rakerLModel)
fitted(rakerSModel)
fitted(rakerCModel) #different values than originally reported because took fitted values from lm not glm, changed doc to glm fitted values
fitted(rakerCModellog)#these are the current values in the doc. But the values in rakerCModel are different because they are not logged
#fitted plots
ggplot(dfraker,aes(x = DOC, y =fitted_L, label= Lake)) + 
  geom_point(aes(colour = Lake)) + 
  geom_smooth(method = 'lm')
ggplot(dfraker,aes(x = DOC, y =fitted_S, label= Lake)) + 
  geom_point(aes(colour = Lake)) + 
  geom_smooth(method = 'lm')
ggplot(dfraker,aes(x = DOC, y =fitted_C, label= Lake)) + 
  geom_point(aes(colour = Lake)) + 
  geom_smooth(method = 'lm')
#Color by basin
ggplot(dfraker,aes(x = DOC, y =fitted_L, label= Lake)) + 
  geom_point(aes(colour = basin.raker)) + 
  geom_smooth(method = 'lm')
ggplot(dfraker,aes(x = DOC, y =fitted_S, label= Lake)) + 
  geom_point(aes(colour = basin.raker)) + 
  geom_smooth(method = 'lm')
ggplot(dfraker,aes(x = DOC, y =fitted_C, label= Lake)) + 
  geom_point(aes(colour = basin.raker)) + 
  geom_smooth(method = 'lm')
-------------------------------------------------------

#Pectoral Fins
setwd("~/")
PecFinDat<-read.csv("PecFinDataNovember.csv")

dfFin<-data.frame(PecFinLength<-PecFinDat$PecFinLengths,
                  lakeID<-PecFinDat$lakeID,
                  fishLengthfin<-PecFinDat$fishStdLength,
                  DOClevel<-PecFinDat$DOCbin,
                  DOC<-PecFinDat$DOC,
                  PecFinBase<-PecFinDat$baseWidth,
                  fishID<-PecFinDat$fishID,
                  finLss<-PecFinDat$finLengthSS,
                  finWss<-PecFinDat$finBaseSS, finLbylakess<-PecFinDat$bylake_finlengthSS,
                  finWbylakess<-PecFinDat$bylake_finwidthSS, finRss<-PecFinDat$ss_ratio, basin.fin<-PecFinDat$basin)
finLModel=lmer(log10(finLss) ~ 1 + log10(DOC) + (1|lakeID), data=dfFin)
finWModel=lmer(log10(finWss) ~ 1 + log10(DOC) + (1|lakeID), data=dfFin)
finRModel=lmer(log10(finRss) ~ 1 + log10(DOC) + (1|lakeID), data=dfFin)
summary(finLModel)
summary(finWModel)
summary(finRModel)
fitted(finLModel)
fitted(finWModel)
fitted(finRModel)


