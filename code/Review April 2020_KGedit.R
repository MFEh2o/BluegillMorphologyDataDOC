library(geomorph)
library(shapes)
library(Morpho)
library(StereoMorph)
library(phangorn)
library(dplyr) # pipes etc
library(lme4) # for mixed models
library(lmerTest)
library(ggplot2) # for plots
library(MuMIn) # for effect sizes
library(here) # for file paths

# Load data ---------------------------------------------------------------
myData_tps <- readland.tps(here("data", "inputs", "FULL_2018_TPS_FILE_UPDATED_09-25-19.TPS"), specID = "imageID")
mylinks <- read.table(here("data", "inputs", "Full_body_links.txt"))
identifiers <- read.table(here("data", "outputs", "Identifiers_Update_2020.txt"), 
                          sep = ",", header = TRUE) # reading in the version that I created from inputs

dimnames(myData_tps)[[3]] <- identifiers$imageID
dimnames(myData_tps)[[1]] <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19") 
GPA.fish <- gpagen(myData_tps, ProcD = TRUE, Proj = TRUE) 
gdf.fish <- geomorph.data.frame(shape=GPA.fish$coords,
                             DOCrange=identifiers$DOCrange,
                             Lake=identifiers$lakeID,
                             DOC=identifiers$lakeDOC,
                             DOCcat=identifiers$DOC,
                             captureMethod=identifiers$captureMethod, 
                             cSize=GPA.fish$Csize, 
                             basin=as.factor(identifiers$Basin))

corePCA <- plotTangentSpace(gdf.fish$shape, groups = as.factor(gdf.fish$DOC), 
                            legend = TRUE)
PCscores <- corePCA$pc.scores[,1:2]
idx <- which(gdf.fish$Lake=="Bay")
Bay <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="Birch")
Birch <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="Crampton")
Crampton <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="Found")
Found <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="Hummingbird")
Hummingbird <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="Little_Crooked")
Little_Crooked <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="Lost")
Lost <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="McCullough")
McCullough <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="Muskellunge")
Muskellunge <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="Oxbow")
Oxbow <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="Papoose")
Papoose <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="Red_Bass")
Red_Bass <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="Squaw")
Squaw <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="Towanda")
Towanda <- PCscores[idx,]


# PC1 plot ----------------------------------------------------------------
plotRefToTarget(corePCA$pc.shapes$PC1min, 
                corePCA$pc.shapes$PC1max,
                method="points", 
                links = mylinks, 
                gridPars = gridPar(tar.pt.bg = "black",
                                   tar.link.col = "black",
                                   tar.link.lwd = 3, 
                                   tar.pt.size = 1, 
                                   pt.size = 1, 
                                   pt.bg = "gray", 
                                   link.lwd = 3), 
                mag=1, 
                useRefPts = TRUE)


# PC2 plot ----------------------------------------------------------------
plotRefToTarget(corePCA$pc.shapes$PC2min, 
                corePCA$pc.shapes$PC2max,
                method="points", 
                links = mylinks, 
                gridPars = gridPar(tar.pt.bg = "black",
                                   tar.link.col = "black",
                                   tar.link.lwd = 3, 
                                   tar.pt.size = 1, 
                                   pt.size = 1, 
                                   pt.bg = "gray", 
                                   link.lwd = 3), 
                mag=1, 
                useRefPts = TRUE)

lm_array <- myData_tps
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
    coor <- coor %*% matrix(c(1,0,0,0,cos(-pi/2),sin(-pi/2), 
                              0,-sin(-pi/2),cos(-pi/2)), 
                            nrow=3, ncol=3)
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

### BACKTRANSFORM MORPHOSPACE ALL INDIVIDUALS
# Create plot box with axes and axis labels
plot(scores[, pcs], type="n", main="Backtransform morphospace",
     xlab=paste0("PC", pcs[1], " (", round(per_var[pcs[1]]), "%)"),
     ylab=paste0("PC", pcs[2], " (", round(per_var[pcs[2]]), "%)"))

# Plot backtransform shapes, changed sign of rotation matrix (resEig$vectors) 
btShapes(scores=scores, vectors=-(resEig$vectors), fcn=plot_fish_lateral,pcs=pcs,
         n=c(4,4), m=dim(lm_array)[2], row.names=dimnames(lm_array)[[1]], 
         pc.margin=c(0.06,0.05), size=0.038, col=gray(0.7))


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

#######################
#BACKTRANSFORM MRPHOSPACE FIGURE WITH MEAN SHAPES PER LAKE
###Mean Shapes Per Lake
lm_array <- myData_tps
gpa_array <- gpagen(myData_tps, ProcD = TRUE, Proj = TRUE)$coords ### takes just the coords

x <- two.d.array(gpa_array) # the following part is necessary to calculate mean shapes per lake
gpa_df = geomorph.data.frame(shape=gpa_array, 
                             Lake= as.factor(identifiers$lakeID))
p <- dim(gpa_array)[1] # the number of landmarks
k <- dim(gpa_array)[2] # the dimensions of the coordinates
Y <- array(NA, dim = c(p,k,length(levels(gpa_df$Lake)))) #new empty array to fill
dimnames(Y)[[3]] <- levels(gpa_df$Lake) # set group levels as new names

for (i in 1:length(levels(gpa_df$Lake))){
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


##RE-RUN STATS FOR SHAPE TO CHECK THEM

#DOC allometric slope
final <- anova(procD.lm(shape ~ log(cSize) + log(DOC) + Lake + log(DOC)*Lake +
                       log(cSize)*log(DOC) ,SS.type = "II", data = gdf.fish),error= c("Residuals","Residuals","Lake","Residuals"))
summary(final)

#Lake allometric slope
final2 <- anova(procD.lm(shape ~ log(cSize) + log(DOC) + Lake + log(DOC)*Lake + 
                       log(cSize)*Lake, SS.type = "II",data = gdf.fish),error= c("Residuals","Residuals","Lake","Residuals"))
summary(final2)

#Slopes dependent on basin and doc
ShapeModelNew <- anova(procD.lm(shape ~ log(cSize) + log(DOC) + basin + Lake + log(cSize):log(DOC) +
                               log(cSize):basin + basin:log(DOC) + log(cSize):basin:log(DOC), SS.type = "II",data = gdf.fish), error= c("Residuals","Residuals","Residuals","Lake", "Residuals", "Residuals", "Residuals"))
summary(ShapeModelNew)

#Shape Model New allows slopes to vary by basin, DOC and basin:DOC 
#I believe this is the closest multivariate shape model to the one we are using for univariate variables. 
#Note, output does not include the interaction between basin and DOC. Not sure why this is ommitted, perhaps because it is explained by three-way interaction?
#Shape Model 2 is exact same output but much shorter to write. 
ShapeModel2 <- anova(procD.lm(shape ~ log(cSize) * basin * log(DOC) + Lake,SS.type = "II",data = gdf.fish),error= c("Residuals","Residuals","Residuals","Lake","Residuals","Residuals","Residuals"))
summary(ShapeModel2)

###### Re-test common slope
fit.DOC.unique<- procD.lm(shape~log(cSize)*log(DOC)*Lake, data = gdf.fish) #unique DOC allometries
fit.DOC.common<- procD.lm(shape~log(cSize)+log(DOC)*Lake, data = gdf.fish) #common DOC allometries
fitcompare <- anova(fit.DOC.common, fit.DOC.unique, error = c("Residuals","Residuals","DOC:Lake"))
summary(fitcompare)
#basin common slope test
fit.basin.unique <- procD.lm(shape~log(cSize)*basin, data = gdf.fish)
fit.basin.common <- procD.lm(shape~log(cSize)+basin, data = gdf.fish)
aov.fit.basin <- anova(fit.basin.common,fit.basin.unique)
summary(aov.fit.basin)

###### Re-test full vs reduced models
reduced_DOC_basin_model <- anova(procD.lm(shape ~ log(cSize) + log(DOC) + basin + Lake + log(DOC)*Lake + log(cSize)*log(DOC) + basin*Lake + log(cSize)*basin,SS.type = "II",data = gdf.fish),error= c("Residuals","Residuals","Residuals","Lake","Residuals","Residuals"))
summary(reduced_DOC_basin_model)  

#Full model below includes size*DOC*basin interaction which is significant so it seems like it is a better choice??
fullmodel <- anova(procD.lm(shape ~ log(cSize) + log(DOC) + basin + Lake + log(DOC)*Lake + log(DOC)*basin + basin*Lake + log(DOC)*Lake*basin + log(cSize)*log(DOC)*basin,SS.type = "II",data = gdf.fish), error=c("Residuals","Residuals","Residuals","Lake","Residuals","Residuals","Residuals"))
summary(fullmodel)
#NOTE:Same result from full model as "SHAPE MODEL NEW" used above and in thesis! So, I think what I did in my thesis was correct

# TPS grids for min and max scores in previous plot
fit <- procD.lm(shape ~ log(cSize) + log(DOC) + basin + Lake  +
                log(cSize)*basin*log(DOC), 
                SS.type = "II", data = gdf.fish)
plot(fit, type = "diagnostics") 

# diagnostic plots, including plotOutliers
plot(fit, type = "diagnostics", outliers = TRUE) 

# PC plot rotated to major axis of fitted values
plot(fit, type = "PC", pch = 19, col = "blue") 
# Regression-type plots

# Use fitted values from the model to make prediction lines
plot(fit, type = "regression", 
     predictor = gdf.fish$cSize, reg.type = "RegScore", 
     pch = 19, col = "green")

# Uses coefficients from the model to find the projected regression scores
fish.plot <- plot(fit, type = "regression", 
                 predictor = gdf.fish$cSize, reg.type = "RegScore", 
                 pch = 21, bg = "yellow") 

preds <- shape.predictor(fit$GM$fitted, x = fish.plot$RegScore, 
                         predmin = min(fish.plot$RegScore), 
                         predmax = max(fish.plot$RegScore))
#To Get Mean Configuration For All (Doesn't really have anything to do with fitted values from model)
M <- GPA.fish$consensus
plotRefToTarget(M, preds$predmin, mag=2)
plotRefToTarget(M, preds$predmax, mag=2)

attributes(fit)
fit$fitted[1:3, ] # the fitted values (first three specimens)
fit$GM$fitted[,, 1:3] # the fitted values as Procrustes coordinates (same specimens)
#Fitted Values Allometry colored by basin
# plotAllometry(fit, size=gdf.fish$cSize, logsz = TRUE, method = "PredLine", pch=19, col=rainbow(2)[gdf.fish$basin], main="Predicted PC1 Values From Model vs. Size (Colored By Basin)") # This was code that was just used for checking--maybe associated with an older version of geomorph? No longer runs, but also not super critical.

#Predicted/Fitted PC1 max and min from Model
plotRefToTarget(preds$predmin, preds$predmax,method="points", links = mylinks, gridPars = gridPar(tar.pt.bg = "black",tar.link.col = "black",tar.link.lwd = 3, tar.pt.size = 1, pt.size = 1, pt.bg = "gray", link.lwd = 3), mag=1, useRefPts = TRUE)

###Univariate Data Check and Run###
#Eyes
dfeye <-read.csv(here("data", "outputs", "eyewidthsFINAL.csv")) # This is the re-created file, generated in recreateEyewidths.R. It does not include by-lake size-standardized values because those values were not used in the subsequent analysis.

# Create a new data frame, changing some of the column names.
dfeye <- dfeye %>%
  select("eyewidth" = eyeWidth, 
         "fishIDeye" = fishID,
         lakeID,
         "fishLength" = fishStdLength,
         DOClevel,
         DOC,
         basin,
         area,
         "maxDepth" = max_Depth,
         captureMethod,
         eyewidth.ss)

# Original Model 
EyeModel <- lmer(log(eyewidth.ss) ~ 1 + log(DOC) + (1|lakeID), data=dfeye)
summary(EyeModel)

# New Model with Basin
EyeModelNew <- lmer(log(eyewidth.ss) ~ 1 + log(DOC) + basin + basin:log(DOC) + (1|lakeID), data = dfeye)
summary(EyeModelNew)

# Get fitted values and add to data frame. Note only one value for each lake, so only 14 different fitted values.
fitted(EyeModel) # fitted values for the first model
fitted(EyeModelNew) # fitted values for the model that includes basins. These are the ones we want to save.

## Add the basin-model fitted values to dfeye
dfeye$fitted <- fitted(EyeModelNew) # one value per lake--only 14 different fitted values.

# Get R2 effect sizes
r.squaredGLMM(EyeModelNew, by_group=TRUE)
#           R2m       R2c
#[1,] 0.02491908 0.3809909 # (ish--may vary slightly)
#R2m is marginal variance which shows variance explained by fixed effects
#R2c is conditional variance which shows variance explained by total model

ggplot(dfeye,aes(x = DOC, y =fitted, label= lakeID)) + 
  geom_point(aes(colour = lakeID))
# Other Plots to Compare (raw eye widths and size standardized eye widths)
ggplot(dfeye,aes(x = DOC, y = eyewidth, label= lakeID)) + 
  geom_point(aes(colour = lakeID))
ggplot(dfeye,aes(x = DOC, y = eyewidth.ss, label= lakeID)) + 
  geom_point(aes(colour = lakeID))

# YESSSSSS all done with the eye widths!!

#Gill Rakers
rakerData <- read.csv("Gill_Rakers_2018_Final.csv")
dfraker<-data.frame(captureMethod<-rakerData$capture_Method,lakeLss=rakerData$avgrakerlengthSS_bylake,
                    lakeSss=rakerData$avgrakerspaceSS_bylake,lakeCss=rakerData$rakercountSS_bylake,
                    rakercount=rakerData$total_RakerNum,DOClevel=rakerData$DOCbin,DOC=rakerData$lakeDOC,
                    Lake=rakerData$lakeID,fishLength=rakerData$fishSL,avgRakerLength=rakerData$avgRakerLength,
                    avgRakerSpace=rakerData$avgRakerSpace,count.ss=rakerData$SS.Count, length.ss=rakerData$SS.Length,
                    space.ss=rakerData$SS.Space,avgL2<-rakerData$avgL_4.7,avgS2<-rakerData$avgS_4.6,
                    avgL2_ss<-rakerData$avgL2_ss,avgS2_ss<-rakerData$avgS2_ss, basin.raker<-rakerData$basin, fitted_L<-rakerData$fitted_L,fitted_S<-rakerData$fitted_S,fitted_C<-rakerData$fitted_C)
#avgL2_ss is size standardizations for average length for rakers 4-7
rakerLModel=lmer(log(avgL2_ss) ~ 1 + log(DOC) + (1|Lake), data=dfraker)
summary(rakerLModel)
rakerSModel=lmer(log(avgS2_ss) ~ 1 + log(DOC) + (1|Lake), data=dfraker)
summary(rakerSModel)
#With Basin
RakerLModelNew=lmer(log(avgL2_ss) ~ 1 + log(DOC) + basin.raker + basin.raker:log(DOC) + (1|Lake), data = dfraker)
summary(RakerLModelNew)
RakerSModelNew=lmer(log(avgS2_ss) ~ 1 + log(DOC) + basin.raker + basin.raker:log(DOC) + (1|Lake), data = dfraker)
summary(RakerSModelNew)
fitted(RakerLModelNew)
fitted(RakerSModelNew)
#use glm for count data
rakerCModel=glmer(log(rakercount) ~ 1 + log(DOC) + (1|Lake),family = poisson,nAGQ=0,control=glmerControl(optimizer = "nloptwrap"), data=dfraker)
summary(rakerCModel)
RakerCModelNew=glmer(rakercount ~ 1 + log(DOC) + basin.raker + basin.raker:log(DOC) + (1|Lake),family = poisson,nAGQ=0,control=glmerControl(optimizer = "nloptwrap"), data=dfraker)
summary(RakerCModelNew)
fitted(RakerCModelNew)
#Still Singular!
#Going to have to use lmer
RakerCModelNew=lmer(rakercount ~ 1 + log(DOC) + basin.raker + basin.raker:log(DOC) + (1|Lake), data = dfraker)
summary(RakerCModelNew)
fitted(RakerCModelNew)
#fitted plots
ggplot(dfraker,aes(x = DOC, y =fitted_L, label= Lake)) + 
  geom_point(aes(colour = Lake)) 
ggplot(dfraker,aes(x = DOC, y =fitted_S, label= Lake)) + 
  geom_point(aes(colour = Lake))
ggplot(dfraker,aes(x = DOC, y =fitted_C, label= Lake)) + 
  geom_point(aes(colour = Lake))

#Other figures to compare
ggplot(dfraker,aes(x = DOC, y =avgL2, label= Lake)) + 
  geom_point(aes(colour = Lake)) 
ggplot(dfraker,aes(x = DOC, y =avgS2, label= Lake)) + 
  geom_point(aes(colour = Lake))
ggplot(dfraker,aes(x = DOC, y =rakercount, label= Lake)) + 
  geom_point(aes(colour = Lake))
ggplot(dfraker,aes(x = DOC, y =avgL2_ss, label= Lake)) + 
  geom_point(aes(colour = Lake)) 
ggplot(dfraker,aes(x = DOC, y =avgS2_ss, label= Lake)) + 
  geom_point(aes(colour = Lake))
ggplot(dfraker,aes(x = DOC, y =count.ss, label= Lake)) + 
  geom_point(aes(colour = Lake))


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
                  finWbylakess<-PecFinDat$bylake_finwidthSS, finRss<-PecFinDat$ss_ratio, 
                  basin.fin<-PecFinDat$basin, fitted_PL=PecFinDat$fitted_PL, fitted_PW=PecFinDat$fitted_PW, fitted_PR=PecFinDat$fitted_PR)
finLModel=lmer(log(finLss)~ 1 + log(DOC) + (1|lakeID), data=dfFin)
# > isSingular(finLModel)
#[1] TRUE
#This means that some "dimensions" of the variance-covariance matrix have been estimated as exactly zero, 
#in output, variance and std. dev. of lake intercepts are both 0
#Getting rid of log infront of finLss gets rid of error. Note,just getting rid of log infront of DOC still results in sing.
finLModel=lmer(finLss~ 1 + log(DOC) + (1|lakeID), data=dfFin)## log removed
finWModel=lmer(log(finWss) ~ 1 + log(DOC) + (1|lakeID), data=dfFin)
finRModel=lmer(log(finRss) ~ 1 + log(DOC) + (1|lakeID), data=dfFin)
summary(finLModel)
summary(finWModel)
summary(finRModel)
fitted(finLModel)
fitted(finWModel)
fitted(finRModel)
#NEW MODELS WITH BASIN
FinLModelNew=lmer(finLss ~ 1 + log(DOC) + basin.fin + basin.fin:log(DOC) + (1|lakeID), data = dfFin)
summary(FinLModelNew) #Still singular had to remove log
#Fin Length now Significant!!
FinWModelNew=lmer(log(finWss) ~ 1 + log(DOC) + basin.fin + basin.fin:log(DOC) + (1|lakeID), data = dfFin)
summary(FinWModelNew)
#Still no effect on base width
FinRModelNew=lmer(log(finRss) ~ 1 + log(DOC) + basin.fin + basin.fin:log(DOC) + (1|lakeID), data = dfFin)
summary(FinRModelNew)
#Still no effect on fin ratio

fitted(FinLModelNew)
fitted(FinWModelNew)
fitted(FinRModelNew)

ggplot(dfFin,aes(x = DOC, y = fitted_PL, label= lakeID)) + 
  geom_point(aes(colour = lakeID)) 

ggplot(dfFin,aes(x = DOC, y = fitted_PW, label= lakeID)) + 
  geom_point(aes(colour = lakeID))

ggplot(dfFin,aes(x = DOC, y = fitted_PR, label= lakeID)) + 
  geom_point(aes(colour = lakeID)) 

ggplot(dfFin,aes(x = DOC, y = PecFinLength, label= lakeID)) + 
  geom_point(aes(colour = lakeID)) 
ggplot(dfFin,aes(x = DOC, y = PecFinBase, label= lakeID)) + 
  geom_point(aes(colour = lakeID))
ggplot(dfFin,aes(x = DOC, y = finLss, label= lakeID)) + 
  geom_point(aes(colour = lakeID)) 
ggplot(dfFin,aes(x = DOC, y = finWss, label= lakeID)) + 
  geom_point(aes(colour = lakeID))
ggplot(dfFin,aes(x = DOC, y = finRss, label= lakeID)) + 
  geom_point(aes(colour = lakeID))

#finAngles
finAngleData<-read.csv("Pec Fin Angles.csv")
dfangle<-data.frame(fishID=finAngleData$fishID,lakeID=finAngleData$lakeID,
                    DOClevel=finAngleData$DOClevel,DOC=finAngleData$DOC,fishLength=finAngleData$fishStdLength,
                    captureMethod<-finAngleData$captureMethod, basin.angle<-finAngleData$basin,
                    finangle<-finAngleData$ang_deg, fittedang<-finAngleData$fitted)

finAModel=lmer(log(finangle) ~ 1 + log(DOC) + (1|lakeID), data=dfangle)
summary(finAModel)

FinAModelNew=lmer(log(finangle) ~ 1 + log(DOC) + basin.angle + basin.angle:log(DOC) + (1|lakeID), data = dfangle)
summary(FinAModelNew)

fitted(FinAModelNew)

ggplot(dfangle,aes(x = DOC, y = fittedang, label= lakeID)) + 
  geom_point(aes(colour = lakeID)) 
#Check to see if relationship with size to see if needs to be size corrected
ggplot(dfangle,aes(x = fishLength, y = finangle, label= lakeID)) + 
  geom_point(aes(colour = lakeID)) 
#No size relationship, can leave it as is.

#Raw Plot
ggplot(dfangle,aes(x = DOC, y = finangle, label= lakeID)) + 
  geom_point(aes(colour = lakeID)) 
