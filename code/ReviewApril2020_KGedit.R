# Main analysis script
# Created by Chelsea Bishop. Edits/reorganization by Kaija Gahm
# Contains main fish size model, as well as models for eye widths, gill raker lengths/spaces, pectoral fin lengths/widths, and pectoral fin insertion angles. Also contains code for figures and tables.


# Load packages -----------------------------------------------------------
# versions can be found in the renv lock file
library(geomorph) # for morphometric analysis
library(shapes) # for morphometric analysis
library(Morpho) # for morphometric analysis
library(StereoMorph) # for morphometric analysis
library(phangorn) # for morphometric analysis?
library(dplyr) # for pipes etc
library(lme4) # for mixed models
library(lmerTest) # to calculate stats on mixed models
library(ggplot2) # for plots
library(MuMIn) # for effect sizes
library(here) # for file paths
library(gridGraphics) # for base R plots
library(ggtext) # for ggplots
source(here("code", "defs.R"))

# Load data ---------------------------------------------------------------
## landmark data
myData_tps <- readland.tps(here("data", "inputs", "FULL_2018_TPS_FILE_UPDATED_09-25-19.TPS"), specID = "imageID")

## manually-defined link data
mylinks <- read.table(here("data", "inputs", "Full_body_links.txt"))

## fish metadata to correspond to the landmarks
identifiers <- read.table(here("data", "outputs", "Identifiers_Update_2020.txt"), 
                          sep = ",", header = TRUE)


# Prepare the data for morphometric analyses -------------------------------
## set imageID's as names for the individual fish (third dimension)
dimnames(myData_tps)[[3]] <- identifiers$imageID

## number the landmarks (first dimension)
dimnames(myData_tps)[[1]] <- as.character(1:19)


# Morphometric analysis ------------------------------------------------------
## perform Procrustes analysis
GPA.fish <- gpagen(myData_tps, ProcD = TRUE, Proj = TRUE) 

## convert results of Procrustes analysis into a usable data frame
gdf.fish <- geomorph.data.frame(shape=GPA.fish$coords,
                                Lake=identifiers$lakeID,
                                DOC=identifiers$lakeDOC,
                                captureMethod=identifiers$captureMethod, 
                                cSize=GPA.fish$Csize, 
                                basin=as.factor(identifiers$Basin))

## PCA scoreplot for overall fish shapes, colored by DOC levels (proxy for lakes)
corePCA <- plotTangentSpace(gdf.fish$shape, 
                            groups = as.factor(gdf.fish$DOC), 
                            legend = TRUE)

## save PC scores to their own object
PCscores <- corePCA$pc.scores[,1:2]

## Split PC scores into individual objects by lake
idx <- which(gdf.fish$Lake=="BA")
Bay <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="BH")
Birch <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="CR")
Crampton <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="FD")
Found <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="HB")
Hummingbird <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="LC")
Little_Crooked <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="LT")
Lost <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="MC")
McCullough <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="MS")
Muskellunge <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="OB")
Oxbow <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="PS")
Papoose <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="RS")
Red_Bass <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="SQ")
Squaw <- PCscores[idx,]

idx <- which(gdf.fish$Lake=="TO")
Towanda <- PCscores[idx,]


# Save plots of the PC1 extremes ------------------------------------------
# Open an pdf file
pdf(here("figures", "fishShapes_pc1_pc2", "pc1Min.pdf"), width = 5, height = 5) 
# plot
plotRefToTarget(corePCA$pc.shapes$PC1min, 
                corePCA$pc.shapes$PC1max,
                method="points", 
                links = mylinks, 
                gridPars = gridPar(tar.pt.bg = "darkorange3",
                                   tar.link.col = "darkorange3",
                                   tar.link.lwd = 3, 
                                   tar.pt.size = 1, 
                                   pt.size = 1, 
                                   pt.bg = "gray", 
                                   link.lwd = 3), 
                mag=1, 
                useRefPts = TRUE)
dev.off() # close the device

# Open an pdf file
pdf(here("figures", "fishShapes_pc1_pc2", "pc1Max.pdf"), width = 5, height = 5) 
# plot
plotRefToTarget(corePCA$pc.shapes$PC1max,
                corePCA$pc.shapes$PC1min, 
                method="points", 
                links = mylinks, 
                gridPars = gridPar(tar.pt.bg = "darkorange3",
                                   tar.link.col = "darkorange3",
                                   tar.link.lwd = 3, 
                                   tar.pt.size = 1, 
                                   pt.size = 1, 
                                   pt.bg = "gray", 
                                   link.lwd = 3), 
                mag=1, 
                useRefPts = TRUE)
dev.off() # close the device

# Save plots of the PC2 extremes ------------------------------------------
# Open an pdf file
pdf(here("figures", "fishShapes_pc1_pc2", "pc2Min.pdf"), width = 5, height = 5) 
# plot
plotRefToTarget(corePCA$pc.shapes$PC2min, 
                corePCA$pc.shapes$PC2max,
                method="points", 
                links = mylinks, 
                gridPars = gridPar(tar.pt.bg = "green3",
                                   tar.link.col = "green3",
                                   tar.link.lwd = 3, 
                                   tar.pt.size = 1, 
                                   pt.size = 1, 
                                   pt.bg = "gray", 
                                   link.lwd = 3), 
                mag=1, 
                useRefPts = TRUE)
dev.off()

# Open an pdf file
pdf(here("figures", "fishShapes_pc1_pc2", "pc2Max.pdf"), width = 5, height = 5) 
# plot
plotRefToTarget(corePCA$pc.shapes$PC2max, 
                corePCA$pc.shapes$PC2min,
                method="points", 
                links = mylinks, 
                gridPars = gridPar(tar.pt.bg = "green3",
                                   tar.link.col = "green3",
                                   tar.link.lwd = 3, 
                                   tar.pt.size = 1, 
                                   pt.size = 1, 
                                   pt.bg = "gray", 
                                   link.lwd = 3), 
                mag=1, 
                useRefPts = TRUE)
dev.off()

# [needs title] ----------------------------------------------------------------
lm_array <- myData_tps
gpa_array <- gpagen(myData_tps, ProcD = TRUE, Proj = TRUE)$coords ### takes just the coords

# Convert array to matrix for PCA
gpa_mat <- t(apply(gpa_array, 3, function(y) matrix(t(y), 1)))

# Perform non-phylogenetic PCA
resEig <- eigen(cov(gpa_mat))

# Get PC scores, changed sign of rotation matrix (resEig$vectors)
scores <- gpa_mat %*% -(resEig$vectors)
#### shapes switched, need to reverse signs to fix

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

# Backtransform morphospace all individuals -------------------------------
# Create plot box with axes and axis labels
plot(scores[, pcs], type="n", main="Backtransform morphospace",
     xlab=paste0("PC", pcs[1], " (", round(per_var[pcs[1]]), "%)"),
     ylab=paste0("PC", pcs[2], " (", round(per_var[pcs[2]]), "%)"))

# Plot backtransform shapes, changed sign of rotation matrix (resEig$vectors) 
btShapes(scores=scores, vectors=-(resEig$vectors), fcn=plot_fish_lateral,pcs=pcs,
         n=c(4,4), m=dim(lm_array)[2], row.names=dimnames(lm_array)[[1]], 
         pc.margin=c(0.06,0.05), size=0.038, col=gray(0.7))

# add points for each lake in a different color
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


# Backtransform morphospace figure with mean shapes per lake ----------
## Mean Shapes Per Lake
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
#### shapes switched need to reverse signs to fix

# Get percent variance explained along each axis
per_var <- (resEig$values / sum(resEig$values))*100

# Make Fig 3 --------------------------------------------------------------
# Open a pdf file
pdf(here("figures", "fishShapes_pc1_pc2", "mainPlot.pdf"), width = 9, height = 9) 
# plot
# Create plot box with axes and axis labels
plot(scores[, pcs], type="n", main="Backtransform morphospace",
     xlab=paste0("PC", pcs[1], " (", round(per_var[pcs[1]]), "%)"),
     ylab=paste0("PC", pcs[2], " (", round(per_var[pcs[2]]), "%)"))


# Plot backtransform shapes, changed sign of rotation matrix (resEig$vectors) 
btShapes(scores=scores, vectors=-(resEig$vectors), 
         fcn=plot_fish_lateral, 
         pcs=pcs, n=c(4,4), m=dim(lm_array)[2], 
         row.names=dimnames(lm_array)[[1]], 
         pc.margin=c(0.06,0.06), 
         size=0.018, 
         col=gray(0.7))
points(scores [,1], scores[,2], col="#41C6EC", pch=19,cex=2)
text(scores [,1], scores [,2], labels=levels(gpa_df$Lake), cex= 0.7)
dev.off()
# For ref-to-Target Figures that go along the axes in figure 3: see pc1Plot and pc2Plot, above.

# Stats for shape ---------------------------------------------------------

#Model structure similar to univariate analyses

#Check that basin is a factor
str(gdf.fish)

#Fit model
#shapeModel <- procD.lm(shape ~ log(cSize) + log(DOC) + basin + log(DOC):basin + Lake, data=gdf.fish, SS.type="II")

#Because this is not a hierarchical model, the log(DOC):basin term is redundant with the Lake term. Drop the former.
shapeModel <- procD.lm(shape ~ log(cSize) + log(DOC) + basin + Lake, data=gdf.fish, SS.type="II",iter=10000)

#Diagnostic plot
plot(shapeModel)

#ANOVA
anova(shapeModel,error=c("Residuals","Residuals","Residuals","Lake"))


# Univariate Data Check and Run -------------------------------------------
# Eye Widths --------------------------------------------------------------
dfeye <-read.csv(here("data", "outputs", "eyewidthsFINAL.csv")) # This is the re-created file, generated in recreateEyewidths.R. It does not include by-lake size-standardized values because those values were not used in the subsequent analysis.

# Create a new data frame, changing some of the column names.
dfeye <- dfeye %>%
  select("eyewidth" = eyeWidth, 
         "fishIDeye" = fishID,
         lakeID,
         "fishLength" = stdLength,
         DOC,
         basin,
         eyewidth.ss)

# Original Model 
# EyeModel <- lmer(log(eyewidth.ss) ~ 1 + log(DOC) + (1|lakeID), data=dfeye)
# summary(EyeModel)

# New Model with Basin
EyeModelNew <- lmer(log(eyewidth.ss) ~ 1 + log(DOC) + basin + basin:log(DOC) + (1|lakeID), data = dfeye)
summary(EyeModelNew)

# Get fitted values and add to data frame. Note only one value for each lake, so only 14 different fitted values.
# fitted(EyeModel) # fitted values for the first model
fitted(EyeModelNew) # fitted values for the model that includes basins. These are the ones we want to save.

# ## Add the basin-model fitted values to dfeye
dfeye$fitted <- fitted(EyeModelNew) # one value per lake--only 14 different fitted values.

# ## Write out a new csv version that includes the fitted values: will be saved in data/outputs. This is so that Chelsea can use the output fitted values in figures etc. if she needs them.
# write.csv(dfeye, file = here("data", "outputs", "eyewidthsFINAL_wFitted.csv"), row.names = F)

# Get R2 effect sizes
r.squaredGLMM(EyeModelNew, by_group = TRUE)
#           R2m       R2c
#[1,] 0.08229317 0.4030357 # (ish--may vary slightly)
# R2m is marginal variance which shows variance explained by fixed effects
# R2c is conditional variance which shows variance explained by total model

# Some plots:
ggplot(dfeye, aes(x = DOC, y = fitted, label = lakeID)) + 
  geom_point(aes(colour = lakeID))
# Other Plots to Compare (raw eye widths and size-standardized eye widths)
ggplot(dfeye,aes(x = DOC, y = eyewidth, label = lakeID)) + 
  geom_point(aes(colour = lakeID))
ggplot(dfeye,aes(x = DOC, y = eyewidth.ss, label = lakeID)) + 
  geom_point(aes(colour = lakeID))

# Gill Rakers -------------------------------------------------------------
dfraker <- read.csv(here("data", "outputs", "Gill_Rakers_2018_Final.csv"))

# Without Basin
# avgL2_ss is size standardizations for average length for rakers 4-7
# rakerLModel <- lmer(log(avgL2_ss) ~ 1 + log(lakeDOC) + (1|lakeID), data=dfraker)
# summary(rakerLModel)
# 
# rakerSModel <- lmer(log(avgS2_ss) ~ 1 + log(lakeDOC) + (1|lakeID), data=dfraker)
# summary(rakerSModel)

# With Basin
## lengths
RakerLModelNew <- lmer(log(avgL2_ss) ~ 1 + log(lakeDOC) + basin + basin:log(lakeDOC) + (1|lakeID), data = dfraker)
summary(RakerLModelNew)

## spaces
RakerSModelNew <- lmer(log(avgS2_ss) ~ 1 + log(lakeDOC) + basin + basin:log(lakeDOC) + (1|lakeID), data = dfraker)
summary(RakerSModelNew)

## save fitted vals
fittedL <- fitted(RakerLModelNew)
fittedS <- fitted(RakerSModelNew)

# GLM for raker count data
## Without Basins
# rakerCModel <- glmer(log(total_RakerNum) ~ 1 + log(lakeDOC) + (1|lakeID),
#                      family = poisson,
#                      nAGQ = 0,
#                      control = glmerControl(optimizer = "nloptwrap"), 
#                      data = dfraker)
# summary(rakerCModel)

## With Basins
# RakerCModelNew = glmer(total_RakerNum ~ 1 + log(lakeDOC) + basin + basin:log(lakeDOC) + (1|lakeID),
#                        family = poisson,
#                        nAGQ = 0,
#                        control = glmerControl(optimizer = "nloptwrap"), 
#                        data = dfraker)
# summary(RakerCModelNew)
# fitted(RakerCModelNew)
#Still Singular! -CB (KG: Well, it isn't anymore, because of the changes I made to the coefficients; see script recreateGillRakers.R for a more detailed explanation.)

# Re-try using lmer instead of glmer
RakerCModelNew <- lmer(total_RakerNum ~ 1 + log(lakeDOC) + 
                         basin + basin:log(lakeDOC) + 
                         (1|lakeID), data = dfraker)
summary(RakerCModelNew) # this is the model we want to use.

## save fitted vals
fittedC <- fitted(RakerCModelNew)

# Add the fitted values to the df ------------------------------------------
dfraker <- dfraker %>%
  mutate(fitted_L = fittedL,
         fitted_S = fittedS,
         fitted_C = fittedC)

# Write out the fitted values to a new sheet ------------------------------
# write.csv(dfraker, file = here("data", "outputs", "Gill_Rakers_2018_Final_wFitted.csv"), row.names = F)

# Make plots of the fitted values -----------------------------------------
# One point per lake

## raker lengths
ggplot(dfraker, aes(x = lakeDOC, y = fitted_L, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) 

## raker spaces
ggplot(dfraker, aes(x = lakeDOC, y = fitted_S, label = lakeID)) +
  geom_point(aes(colour = lakeID))

## raker counts
ggplot(dfraker, aes(x = lakeDOC, y = fitted_C, label = lakeID)) +
  geom_point(aes(colour = lakeID))

# Other figures to compare ------------------------------------------------
# Original values
## lengths
ggplot(dfraker, aes(x = lakeDOC, y = avgL_4.7, label= lakeID)) + 
  geom_point(aes(colour = lakeID)) 

## spaces
ggplot(dfraker, aes(x = lakeDOC, y = avgS_4.6, label= lakeID)) + 
  geom_point(aes(colour = lakeID))

## counts
ggplot(dfraker, aes(x = lakeDOC, y = total_RakerNum, label= lakeID)) + 
  geom_point(aes(colour = lakeID))

# Size-standardized values
## lengths
ggplot(dfraker, aes(x = lakeDOC, y = avgL2_ss, label= lakeID)) + 
  geom_point(aes(colour = lakeID)) 

## spaces
ggplot(dfraker, aes(x = lakeDOC, y = avgS2_ss, label= lakeID)) + 
  geom_point(aes(colour = lakeID))

## size-standardization is not applicable for counts.

# Pectoral Fins -----------------------------------------------------------
dfFin <- read.csv(here("data", "outputs", "PecFinDataNovemberFINAL.csv"))

## Add basins
nrow(dfFin) # 218
dfFin <- dfFin %>%
  left_join(lakeInfo %>%
              select(lakeID, basin),
            by = "lakeID")
nrow(dfFin) # 218

finLModel <- lmer(log(finLengthSS)~ 1 + log(DOC) + (1|lakeID), data = dfFin) # XXX Chelsea: since this model is no longer singular, can we keep the log?

finLModel <- lmer(finLengthSS ~ 1 + log(DOC) + (1|lakeID), data = dfFin) ## log removed
finWModel <- lmer(log(finBaseSS) ~ 1 + log(DOC) + (1|lakeID), data = dfFin)
finRModel <- lmer(log(finRatioSS) ~ 1 + log(DOC) + (1|lakeID), data = dfFin)

summary(finLModel)
summary(finWModel)
summary(finRModel)

fitted(finLModel)
fitted(finWModel)
fitted(finRModel)

# New models with basin
FinLModelNew <- lmer(finLengthSS ~ 1 + log(DOC) + basin + basin:log(DOC) + 
                    (1|lakeID), data = dfFin)
summary(FinLModelNew) #Still singular had to remove log # XXX the log version is no longer singular--can we add the log back in?

FinWModelNew <- lmer(log(finBaseSS) ~ 1 + log(DOC) + basin + basin:log(DOC) + 
                       (1|lakeID), data = dfFin)
summary(FinWModelNew)
#Still no effect on base width  XXX does this comment still make sense?

FinRModelNew <- lmer(log(finRatioSS) ~ 1 + log(DOC) + basin + 
                       basin:log(DOC) + (1|lakeID), data = dfFin)
summary(FinRModelNew)
#Still no effect on fin ratio # XXX does this comment still make sense?

dfFin$fitted_PL <- fitted(FinLModelNew)
dfFin$fitted_PW <- fitted(FinWModelNew)
dfFin$fitted_PR <- fitted(FinRModelNew)

# Plots of fitted values
ggplot(dfFin, aes(x = DOC, y = fitted_PL, label= lakeID)) + 
  geom_point(aes(colour = lakeID)) 

ggplot(dfFin, aes(x = DOC, y = fitted_PW, label= lakeID)) + 
  geom_point(aes(colour = lakeID))

ggplot(dfFin, aes(x = DOC, y = fitted_PR, label= lakeID)) + 
  geom_point(aes(colour = lakeID)) 

# Plots of non-fitted values
## raw
ggplot(dfFin, aes(x = DOC, y = pecFinLength, label= lakeID)) + 
  geom_point(aes(colour = lakeID)) 

ggplot(dfFin,aes(x = DOC, y = pecFinBaseWidth, label= lakeID)) + 
  geom_point(aes(colour = lakeID))

## size-standardized
ggplot(dfFin,aes(x = DOC, y = finLengthSS, label= lakeID)) + 
  geom_point(aes(colour = lakeID)) 

ggplot(dfFin,aes(x = DOC, y = finBaseSS, label= lakeID)) + 
  geom_point(aes(colour = lakeID))

ggplot(dfFin,aes(x = DOC, y = finRatioSS, label= lakeID)) + 
  geom_point(aes(colour = lakeID))

# Pec Fin Angles ----------------------------------------------------------
dfangle <- read.csv(here("data", "outputs", "PecFinAnglesFINAL.csv")) %>%
  rename("finangle" = "pecFinInsertionAngle")

# Without basin
# finAModel <- lmer(log(finangle) ~ 1 + log(DOC) + (1|lakeID), data = dfangle)
# summary(finAModel)

# With basin
FinAModelNew <- lmer(log(finangle) ~ 1 + log(DOC) + basin + basin:log(DOC) + (1|lakeID),
                     data = dfangle)
summary(FinAModelNew)

# Get fitted values and add to data frame. 
fitted(FinAModelNew) # yep, these are the right values to use.

dfangle$fitted <- fitted(FinAModelNew)

## Write out a new csv version that includes the fitted values: will be saved in data/outputs. This is so that Chelsea can use the output fitted values in figures etc. if she needs them.
# write.csv(dfangle, file = here("data", "outputs", "PecFinAnglesFINAL_wFitted.csv"), row.names = F)

# Plots:
ggplot(dfangle, aes(x = DOC, y = fitted, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) 

# Raw Plot
ggplot(dfangle, aes(x = DOC, y = finangle, label= lakeID)) + 
  geom_point(aes(colour = lakeID)) 

# Figures -----------------------------------------------------------------
## Code reorganized and expanded by Kaija

# Define manual vector of lake shapes to use in the plots
# we want shape 21 (circle) for great lakes (4) and shape 22 (square) for mississippi (7)
lakeShapes <- lakeInfo %>%
  select(lakeID, DOC, basin) %>%
  arrange(-DOC) %>%
  mutate(shape = case_when(basin == 4 ~ 21,
                           basin == 7 ~ 22)) %>%
  pull(shape)

# Define a few colors to create a gradient through. We have to make sure that brown maps to higher DOC values and light blue maps to lower DOC values.

# Figures 1, 2, 3 ---------------------------------------------------------
## Fig 1.: Map
### I don't think we need code for this one. Can just use the same image she had before.

## Fig 2: Landmarks
### Likewise, no code needed; can use same image.

## Fig 3. PC1 vs PC2 plot on DOC gradient, fish body shape. 
### Saved components of this figure above, in sections "PC1 plots" and "PC2 plots"

# Figure 4 ----------------------------------------------------------------
## Pec fins vs. DOC, with gradient. Panels: A (pec fin length), B (pec fin base width), C (length:width ratio), C (insertion anlge)

# Figure 5 ----------------------------------------------------------------
## Gill rakers vs. DOC, with gradient. Panels: A (average gill raker length), B (average gill raker space), C (Gill raker total number)
## some modifications to the data so it plots properly
## Fig. 5. Gill rakers vs. DOC, with gradient. Panels: A (average gill raker length), B (average gill raker space), C (Gill raker total number)
dfraker <- dfraker %>%
  mutate(lakeID = stringr::str_replace(lakeID, "_", " "),
         basin = case_when(basin == "4" ~ "Great Lakes Watershed",
                           basin == "7" ~ "Mississippi Watershed"),
         lakeID = fct_reorder(lakeID, lakeDOC, .desc = T))

## Exponentiate the fitted values to convert them back to real mm units # 
## make panel A, with legend (we'll remove the legend when we plot it)
panelA <- dfraker %>%
  ggplot(aes(x = lakeDOC, y = exp(fitted_L)))+
  geom_point(aes(fill = lakeID, shape = basin), size = 5, 
             col = "black", stroke = 1)+
  theme_classic()+
  scale_fill_manual(name = "**Lakes**",
                    values = lakeColors,
                    guide = guide_legend(override.aes = list(shape = lakeShapes))) +
  scale_shape_manual(name = "",
                     values = c(21, 22))+
  labs(y = "Average Gill Raker Length (mm)")+
  theme(legend.title = element_markdown(),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_blank())

## make panel B
panelB <- dfraker %>%
  ggplot(aes(x = lakeDOC, y = exp(fitted_S)))+
  geom_point(aes(fill = lakeID, shape = basin), size = 5, 
             col = "black", stroke = 1)+
  theme_classic()+
  scale_fill_manual(name = "**Lakes**",
                    values = lakeColors,
                    guide = guide_legend(override.aes = list(shape = lakeShapes))) +
  scale_shape_manual(name = "",
                     values = c(21, 22))+
  labs(y = "Average Gill Raker Space Width (mm)")+
  theme(legend.title = element_markdown(),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_blank(),
        legend.position = "none")

## Make panel C, with x-axis label
### not exponentiating the fitted values on this one because they weren't log-transformed.
panelC <- dfraker %>%
  ggplot(aes(x = lakeDOC, y = fitted_C))+
  geom_point(aes(fill = lakeID, shape = basin), size = 5, 
             col = "black", stroke = 1)+
  theme_classic()+
  scale_fill_manual(name = "**Lakes**",
                    values = lakeColors,
                    guide = guide_legend(override.aes = list(shape = lakeShapes))) +
  scale_shape_manual(name = "",
                     values = c(21, 22))+
  labs(x = "Dissolved Organic Carbon (mg/L)",
       y = "Gill Raker number")+
  theme(legend.title = element_markdown(),
        text = element_text(family = "Helvetica"),
        legend.position = "none")

## extract the legend from panel A
legend <- cowplot::get_legend(panelA)

## make the three stacked plots
allPanels <- cowplot::plot_grid(panelA + 
                                 theme(legend.position = "none"),
                               panelB, panelC, ncol = 1, align = "v")

## Save the main plot panels
pdf(here("figures", "gillRakers", "gillRakersPlotPanels.pdf"), height = 13, width = 5)
allPanels
dev.off()

## Save the legend separately
legendPlot <- cowplot::plot_grid(legend)
pdf(here("figures", "gillRakers", "legend.pdf"), height = 10, width = 4)
legendPlot
dev.off()

# Make the gradient legend alone
gradientLegend <- cowplot::plot_grid(cowplot::get_legend(dfraker %>%
  ggplot(aes(x = lakeDOC, y = avgS_4.6, col = lakeDOC))+
  geom_point()+
  scale_color_gradientn(colors = lakeColors, name = "DOC")))

## Save the gradient legend alone in the main figures folder
pdf(here("figures", "gradientLegend.pdf"), height = 10, width = 4)
gradientLegend
dev.off()
#XXX START HERE

# Figure 6 ----------------------------------------------------------------
# Eye width across DOC gradient. One panel.
## exponentiating these values to undo the log-transformation
eyePlot <- dfeye %>%
  mutate(lakeID = stringr::str_replace(lakeID, "_", " ")) %>%
  mutate(basin = case_when(basin == "4" ~ "Great Lakes Watershed",
                           basin == "7" ~ "Mississippi Watershed"),
         lakeID = factor(lakeID, levels = names(lakeColors))) %>%
  ggplot(aes(x = DOC, y = exp(fitted)))+
  geom_point(aes(fill = lakeID, shape = basin), size = 5, 
             col = "black", stroke = 1)+
  theme_classic()+
  scale_fill_manual(name = "**Lakes**",
                    values = lakeColors,
                    guide = guide_legend(override.aes = list(shape = lakeShapes))) +
  scale_shape_manual(name = "",
                     values = c(21, 22))+
  labs(x = "Dissolved Organic Carbon (mg/L)",
       y = "Eye Width (mm)")+
  theme(legend.title = element_markdown(),
        text = element_text(family = "Helvetica"),
        legend.position = "none")

## Save the main plot panel
pdf(here("figures", "eyeWidths", "eyeWidthsPlotPanel.pdf"), height = 4, width = 7)
eyePlot
dev.off()

## Save the legend separately (same legend as in the previous plot, scaled down a bit)
pdf(here("figures", "eyeWidths", "legend.pdf"), height = 4, width = 2)
legendPlot
dev.off()

# Tables ------------------------------------------------------------------

# Table 1 -----------------------------------------------------------------
## Summary of lake characteristics and sampling for survey lakes. DOC is the mean dissolved organic carbon concentration.
### Don't need to re-create this table: it should be the same.

# Table 2 -----------------------------------------------------------------
## Multivariate analysis of covariance for bluegill shape data. Significant results are bold.
