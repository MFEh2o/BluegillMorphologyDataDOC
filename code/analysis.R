# Final analysis code for the Bluegill Morphology project
# Created by Kaija Gahm in May 2021
# Adapted from ReviewApril2020_KGedit.R--this is a cleaned-up version, with no substantive changes to which analyses are performed.
# Contains main fish size model, as well as models for eye widths, gill raker lengths/spaces, pectoral fin lengths/widths, and pectoral fin insertion angles. Also contains code for figures and tables.

# Load packages -----------------------------------------------------------
# versions can be found in the renv lock file
library(geomorph) # for morphometric analysis (used)
library(StereoMorph) # for morphometric analysis (used)
library(dplyr) # for pipes etc
library(lme4) # for mixed models
library(ggplot2) # for plots
library(MuMIn) # for model effect sizes
library(here) # for file paths
library(gridGraphics) # for base R plots
library(ggtext) # for markdown titles and captions
library(cowplot) # for putting figures together at the end
source(here("code", "defs.R")) # additional variable/function definitions
source(here("code", "plotFishLateral.R"))

# Load data ---------------------------------------------------------------
## Landmark data
landmarks <- readland.tps("kaijagahm/Desktop/BluegillMorphologyDataDOC/data/inputs/FULL_2018_TPS_FILE_UPDATED_09-25-19.TPS", 
                          specID = "imageID")

landmarks <- readland.tps(here("data", "inputs", 
                               "FULL_2018_TPS_FILE_UPDATED_09-25-19.TPS"), 
                          specID = "imageID")

## Manually-defined link data
links <- read.table(here("data", "inputs", "Full_body_links.txt"))

## Fish metadata to correspond to the landmarks
identifiers <- read.table(here("data", "outputs", 
                               "Identifiers_Update_2020.txt"), 
                          sep = ",", header = TRUE)

# Landmark repeatability analysis -----------------------------------------
# Following the procedure outlined on page 258 of Zelditch 2012, and in this blog post: https://www.r-bloggers.com/2015/04/tips-tricks-8-examining-replicate-error/

# 1. Read the coordinate data into R (and format it)
landmarks_replicates <- readland.tps(here("data", "inputs", 
                                          "Total Replicates Coords.TPS"),
                                     specID = "imageID")

### Make data frame
gdf.rep <- geomorph.data.frame(shape = landmarks_replicates)

# 2. Use gpagen() to perform a Procrustes Superimposition
replicates_superimp <- gpagen(landmarks_replicates, ProcD = TRUE, Proj = TRUE)  
plot(replicates_superimp) 
# all of these look reasonable--no landmarks way far away from where they should be.
replicates_corePCA <- plotTangentSpace(replicates_superimp$coords)

### Classifiers
classifiers <- read.table(here("data", "outputs", "Replicates Identifiers.txt"), 
                          sep = ",", header = TRUE)
gdf.rep$replic <- as.factor(classifiers$replicateNo)
gdf.rep <- geomorph.data.frame(shape = landmarks_replicates, 
                               replic = classifiers$replicateNo, 
                               fishID = classifiers$imageID,
                               cSize = replicates_superimp$Csize)

# 3. Perform a Procrustes ANOVA
### ANOVA
replicates_anova <- procD.lm(shape ~ fishID + fishID:replic, data = gdf.rep) 
summary(replicates_anova)

###             Df     SS      MS     Rsq        F       Z Pr(>SS)   
#fishID         59 602102 10205.1 0.97561 187.1758 10.3693   0.001 **
#fishID:replic  60   8509   141.8 0.01379   2.6012  5.3186   0.001 **
#Residuals     120   6543    54.5 0.01060                            
#Total         239 617154        

# To calculate the repeatability of our digitizing ability, we subtract the MS of the replicate term from the individual term and divide by 4 (because we measured each fish 4 times):
#  (MS(fishID) – MS(fishID:replic))/4
# To make these calculations easier, I'm going to extract the MS terms of the ANOVA output as their own variables with descriptive names: 
ms_inter <- replicates_anova$aov.table[1,3] # MS term for the among-fish variability
ms_intra <- replicates_anova$aov.table[2,3] # MS term for the within-individual, among-replicates variability
part1 <- (ms_inter - ms_intra)/4

# Then we calculate the ratio of this value (part1) to the total MS:
# ((MS(fishID) – MS(fishID:replic))/4 ) / (MS(fishID) + MS(fishID:replic)) 

repeatability <- part1/(ms_intra + part1)
# repeatability is 0.947.

# Prepare non-replicate data for morphometric analyses ----------------------
## Set imageID's as names for the individual fish (third dimension)
dimnames(landmarks)[[3]] <- identifiers$imageID

## Number the landmarks (first dimension)
dimnames(landmarks)[[1]] <- as.character(1:19)

# Morphometric analysis ------------------------------------------------------
## Perform Procrustes analysis
gpa <- gpagen(landmarks, ProcD = TRUE, Proj = TRUE) 

## convert results of Procrustes analysis into a usable data frame
gdf <- geomorph.data.frame(shape = gpa$coords,
                           Lake = identifiers$lakeID,
                           DOC = identifiers$lakeDOC,
                           captureMethod = identifiers$captureMethod, 
                           cSize = gpa$Csize, 
                           basin = as.factor(identifiers$Basin))

## PCA scoreplot for overall fish shapes, colored by DOC levels (proxy for lakes)
fish_corePCA <- plotTangentSpace(gdf$shape, 
                                 groups = as.factor(gdf$DOC), 
                                 legend = TRUE)

## save PC scores to their own object
fish_pcScores <- fish_corePCA$pc.scores[,1:2]

## Split PC scores into individual list elements by lake
lakePCScores <- lapply(lakesHighLow, function(x){
  fish_pcScores[which(gdf$Lake == x),]
}) %>%
  setNames(., lakesHighLow)

# Take just the coords from the gpa object
gpa_coords <- gpa$coords 

# Convert array to matrix for PCA
gpa_mat_indiv <- t(apply(gpa_coords, 3, function(y) matrix(t(y), 1)))

# Perform non-phylogenetic PCA
resEig_indiv <- eigen(cov(gpa_mat_indiv))

# Get PC scores, changed sign of rotation matrix (resEig_indiv$vectors)
scores_indiv <- gpa_mat_indiv %*% -(resEig_indiv$vectors)
#### shapes switched, need to reverse signs to fix

# Get percent variance explained along each axis
per_var_indiv <- (resEig_indiv$values / sum(resEig_indiv$values))*100

# We're going to plot PC's 1 and 2
pcs <- 1:2

# (all individuals) Backtransform morphospace all individuals -------------------------------
# I've created three versions of this plot, all on the axes defined from a PCA on individual points (not lake means). One plot shows the lake means unlabeled, another shows the labeled lake means, and the last shows all individual points.

# Create data frame of means by lake, for use later in plotting the legend
means <- lapply(lakePCScores, function(x) data.frame(mn1 = mean(x[,1]),
                                                     mn2 = mean(x[,2]))) %>%
  data.table::rbindlist() %>%
  as.data.frame()
means$lakeID <- names(lakePCScores)

# 1. PLOT WITH MEAN POINTS ONLY, no legend
# Open a pdf file
pdf(here("figures", "fishShapes_pc1_pc2", "mainPlot_meansOnly_nolegend.pdf"), width = 9.5, height = 7) 
# Create plot box with axes and axis labels
par(mar = c(6, 6, 6, 6)) # set margins around the plot so we can read the axis labels
plot(scores_indiv[, pcs], type = "n", # type = "n" plots the axes without plotting the scores, so we can add the scores later
     main = "Backtransform morphospace", # title
     xlab = paste0("PC", pcs[1], " (", round(per_var_indiv[pcs[1]]), "%)"),
     ylab = paste0("PC", pcs[2], " (", round(per_var_indiv[pcs[2]]), "%)"),
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5) # the cex. commands scale up the text to make it larger and more easily readable

# Plot backtransform shapes, changed sign of rotation matrix (resEig_indiv$vectors) 
# There's no documentation for the btShapes function (i.e. ?btShapes or help(btShapes) won't get you anything), but see this tutorial (https://aaronolsen.github.io/tutorials/morphometrics/backtransform.html) for an explanation of the input parameters to btShapes. 
# btShapes_wrapper function is defined in defs.R
btShapes_wrapper(sc = scores_indiv, vc = -(resEig_indiv$vectors))

# add mean points for each lake in a different color
for(i in 1:nrow(means)){
  points(means[i, 1], means[i, 2], 
         bg = lakeColorsHighLow[i], 
         col = "black",
         lwd = 2,
         pch = lakeShapesHighLow[i], cex = 3)
}
dev.off()

# 2. PLOT WITH MEAN POINTS ONLY, with legend
# Open a pdf file
pdf(here("figures", "fishShapes_pc1_pc2", "mainPlot_meansOnly.pdf"), width = 9.5, height = 7) 
# Create plot box with axes and axis labels
par(mar = c(6, 6, 6, 6)) # set margins around the plot so we can read the axis labels
plot(scores_indiv[, pcs], type = "n", # type = "n" plots the axes without plotting the scores, so we can add the scores later
     main = "Backtransform morphospace", # title
     xlab = paste0("PC", pcs[1], " (", round(per_var_indiv[pcs[1]]), "%)"),
     ylab = paste0("PC", pcs[2], " (", round(per_var_indiv[pcs[2]]), "%)"),
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5) # the cex. commands scale up the text to make it larger and more easily readable

btShapes_wrapper(sc = scores_indiv, vc = -(resEig_indiv$vectors))

# add mean points for each lake in a different color
for(i in 1:nrow(means)){
  points(means[i, 1], means[i, 2], 
         bg = lakeColorsHighLow[i], 
         col = "black",
         lwd = 2,
         pch = lakeShapesHighLow[i], cex = 3)
}

legend("topright", legend = lakesHighLow, 
       pch = 19, col = lakeColorsHighLow)
dev.off()

# (per-lake means) Backtransform morphospace figure with mean shapes per lake ----------
## Mean Shapes Per Lake
x <- two.d.array(gpa_coords) # the following part is necessary to calculate mean shapes per lake
gdf_lakes <- geomorph.data.frame(shape = gpa_coords, 
                                 Lake = as.factor(identifiers$lakeID))
p <- dim(gpa_coords)[1] # the number of landmarks
k <- dim(gpa_coords)[2] # the dimensions of the coordinates
gpa_coords_lakes <- array(NA, dim = c(p, k, length(levels(gdf_lakes$Lake)))) #new empty array to fill
dimnames(gpa_coords_lakes)[[3]] <- levels(gdf_lakes$Lake) # set group levels as new names

for (i in 1:length(levels(gdf_lakes$Lake))){
  grp <- x[which(gdf_lakes$Lake == levels(gdf_lakes$Lake)[i]),]
  foo <- arrayspecs(grp, p, k)
  gpa_coords_lakes[,,i] <- mshape(foo) # place into the new 3D array
}

# Convert array to matrix for PCA
gpa_mat_lakes <- t(apply(gpa_coords_lakes, 3, function(y) matrix(t(y), 1)))

# Perform non-phylogenetic PCA
resEig_lakes <- eigen(cov(gpa_mat_lakes))

# Get PC scores, changed sign of rotation matrix (resEig$vectors)
scores_lakes <- gpa_mat_lakes %*% -(resEig_lakes$vectors)

# Get percent variance explained along each axis
per_var_lakes <- (resEig_lakes$values / sum(resEig_lakes$values))*100

# Plot
# Create plot box with axes and axis labels
par(mar = c(6, 6, 6, 6)) 
plot(scores_lakes[, pcs], type = "n", main = "Backtransform morphospace",
     xlab = paste0("PC", pcs[1], " (", round(per_var_lakes[pcs[1]]), "%)"),
     ylab = paste0("PC", pcs[2], " (", round(per_var_lakes[pcs[2]]), "%)"), 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)


# Plot backtransform shapes, changed sign of rotation matrix (resEig$vectors) 
btShapes_wrapper(sc = scores_lakes, vc = -(resEig_lakes$vectors),
                 sz = 0.018)

points(scores_lakes[,1], scores_lakes[,2], col = "#41C6EC", pch = 19, cex = 2)
text(scores_lakes[,1], scores_lakes[,2], labels = levels(gdf_lakes$Lake), cex = 1.5)
# For ref-to-Target Figures that go along the axes in figure 3: see pc1Plot and pc2Plot, above.

# Save plots of the PC1 extremes ------------------------------------------
# Open a pdf file
pdf(here("figures", "fishShapes_pc1_pc2", "pc1Min.pdf"), width = 5, height = 5) 
# plot
plotRefToTarget(
  # reference = max (gray)
  fish_corePCA$pc.shapes$PC1max, 
  # target = min (orange)
  fish_corePCA$pc.shapes$PC1min, 
  method = "points", 
  links = links, 
  gridPars = gridPar(tar.pt.bg = "darkorange3",
                     tar.link.col = "darkorange3",
                     tar.link.lwd = 3, 
                     tar.pt.size = 1, 
                     pt.size = 1, 
                     pt.bg = "gray", 
                     link.lwd = 3), 
  mag = 1, 
  useRefPts = TRUE)
dev.off() # close the device

# Open a pdf file
pdf(here("figures", "fishShapes_pc1_pc2", "pc1Max.pdf"), width = 5, height = 5) 
# plot
plotRefToTarget(
  # reference = min (gray)
  fish_corePCA$pc.shapes$PC1min,
  # target = max (orange)
  fish_corePCA$pc.shapes$PC1max,
  method="points", 
  links = links, 
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
# Open a pdf file
pdf(here("figures", "fishShapes_pc1_pc2", "pc2Min.pdf"), width = 5, height = 5) 
# plot
plotRefToTarget(
  # reference = max (gray)
  fish_corePCA$pc.shapes$PC2max,
  # target = min (green)
  fish_corePCA$pc.shapes$PC2min,
  method="points", 
  links = links, 
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

# Open a pdf file
pdf(here("figures", "fishShapes_pc1_pc2", "pc2Max.pdf"), width = 5, height = 5) 
# plot
plotRefToTarget(
  # reference = min (gray)
  fish_corePCA$pc.shapes$PC2min,
  # target = max (green)
  fish_corePCA$pc.shapes$PC2max,
  method="points", 
  links = links, 
  gridPars = gridPar(tar.pt.bg = "green3",
                     tar.link.col = "green3",
                     tar.link.lwd = 3, 
                     tar.pt.size = 1, 
                     pt.size = 1, 
                     pt.bg = "gray", 
                     link.lwd = 3), 
  mag = 1, 
  useRefPts = TRUE)
dev.off()

# Stats for shape ---------------------------------------------------------
# Model structure similar to univariate analyses
# Check that basin is a factor
class(gdf$basin)

# Fit model
# shapeModel <- procD.lm(shape ~ log(cSize) + log(DOC) + basin + log(DOC):basin + Lake, data = gdf, SS.type = "II")

# Because this is not a hierarchical model, the log(DOC):basin term is redundant with the Lake term. Drop the former.
shapeModel <- procD.lm(shape ~ log(cSize) + log(DOC) + basin + Lake, 
                       data = gdf, SS.type = "II", iter = 10000)

# Diagnostic plot
plot(shapeModel)

#ANOVA
anova(shapeModel, error = c("Residuals", "Residuals", "Residuals", "Lake"))

# Univariate Data Check and Run -------------------------------------------
# Eye Widths --------------------------------------------------------------
dfeye <- read.csv(here("data", "outputs", "eyewidthsFINAL.csv")) 

# Create a new data frame, changing some of the column names.
dfeye <- dfeye %>%
  select("eyewidth" = eyeWidth, 
         "fishIDeye" = fishID,
         lakeID,
         "fishLength" = stdLength,
         DOC,
         basin,
         eyewidth.ss)

# Check structure of dfeye
str(dfeye)

# Set basin to factor
dfeye$basin <- as.factor(dfeye$basin)

# New Model with Basin
EyeModelNew <- lmer(log(eyewidth.ss) ~ 1 + log(DOC) + basin + 
                      basin:log(DOC) + (1|lakeID), data = dfeye)
summary(EyeModelNew)

# Calculate profile CIs for parameter estimates
confEyeModelNew <- confint(EyeModelNew, level = 0.95, method = 'profile')

# Assemble parameter estimates, CIs, R2
tableEyeModelNew <- makeSummary(model = EyeModelNew, 
                                conf = confEyeModelNew)
tableEyeModelNew

# Add the fitted values to dfeye
dfeye$fitted <- fitted(EyeModelNew) 

# Some plots:
ggplot(dfeye, aes(x = DOC, y = fitted, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) # fitted values
ggplot(dfeye,aes(x = DOC, y = eyewidth, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) # raw eyeWidth values
ggplot(dfeye,aes(x = DOC, y = eyewidth.ss, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) # size-standardized eyeWidth values

# Gill Rakers -------------------------------------------------------------
dfraker <- read.csv(here("data", "outputs", "Gill_Rakers_2018_Final.csv"))

# Check structure of dfraker
str(dfraker)

# Set basin to factor
dfraker$basin <- as.factor(dfraker$basin)

# avgL2_ss is size standardization for average length of rakers 4-7
## lengths
RakerLModelNew <- lmer(log(avgL2_ss) ~ 1 + log(lakeDOC) + basin + 
                         basin:log(lakeDOC) + (1|lakeID), data = dfraker)
summary(RakerLModelNew)

# avgS2_ss is size standardization for average space between rakers, for rakers 4-6
## spaces
RakerSModelNew <- lmer(log(avgS2_ss) ~ 1 + log(lakeDOC) + basin + 
                         basin:log(lakeDOC) + (1|lakeID), data = dfraker)
summary(RakerSModelNew)

# Calculate profile CIs for parameter estimates
confRakerLModelNew <- confint(RakerLModelNew, level = 0.95, method = 'profile')
confRakerSModelNew <- confint(RakerSModelNew, level = 0.95, method = 'profile')

# Assemble parameter estimates, CIs, R2
tableRakerLModelNew <- makeSummary(model = RakerLModelNew, 
                                   conf = confRakerLModelNew)
tableRakerSModelNew <- makeSummary(model = RakerSModelNew, 
                                   conf = confRakerSModelNew)

## save fitted vals
fittedL <- fitted(RakerLModelNew)
fittedS <- fitted(RakerSModelNew)

# GLM for raker count data
RakerCModelNew <- lmer(log(total_RakerNum) ~ 1 + log(lakeDOC) + basin + 
                         basin:log(lakeDOC) + (1|lakeID), data = dfraker)
summary(RakerCModelNew)

# Calculate profile CIs for parameter estimates
confRakerCModelNew <- confint(RakerCModelNew, level = 0.95, method = 'profile')

# Assemble parameter estimates, CIs, R2
tableRakerCModelNew <- makeSummary(model = RakerCModelNew, 
                                   conf = confRakerCModelNew) 

## save fitted vals
fittedC <- fitted(RakerCModelNew)

# Add the fitted values to the df ------------------------------------------
dfraker <- dfraker %>%
  mutate(fitted_L = fittedL,
         fitted_S = fittedS,
         fitted_C = fittedC)

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

## (size-standardization is not applicable for counts.)

# Pectoral Fins -----------------------------------------------------------
dfFin <- read.csv(here("data", "outputs", "PecFinDataNovemberFINAL.csv"))

## Add basins
nrow(dfFin) # 218
dfFin <- dfFin %>%
  left_join(lakeInfo %>%
              select(lakeID, basin),
            by = "lakeID")
nrow(dfFin) # still 218, good

# Check structure of dfFin
str(dfFin)

# Set basin to factor
dfFin$basin <- as.factor(dfFin$basin)

# New models with basin
FinLModelNew <- lmer(log(finLengthSS) ~ 1 + log(DOC) + basin + 
                       basin:log(DOC) + (1|lakeID), data = dfFin)
summary(FinLModelNew)
plot(fitted(FinLModelNew)~dfFin$DOC,pch=as.numeric(dfFin$basin))

FinWModelNew <- lmer(log(finBaseSS) ~ 1 + log(DOC) + basin + 
                       basin:log(DOC) + (1|lakeID), data = dfFin)
summary(FinWModelNew)

FinRModelNew <- lmer(log(finRatioSS) ~ 1 + log(DOC) + basin + 
                       basin:log(DOC) + (1|lakeID), data = dfFin)
summary(FinRModelNew)

# Calculate profile CIs for parameter estimates
confFinLModelNew <- confint(FinLModelNew, level = 0.95, method = 'profile')
confFinWModelNew <- confint(FinWModelNew, level = 0.95, method = 'profile')
confFinRModelNew <- confint(FinRModelNew, level = 0.95, method = 'profile')

# Assemble parameter estimates, CIs, R2
tableFinLModelNew <- makeSummary(model = FinLModelNew, 
                                 conf = confFinWModelNew)
tableFinWModelNew <- makeSummary(model = FinWModelNew, 
                                 conf = confFinWModelNew)
tableFinRModelNew <- makeSummary(model = FinRModelNew, 
                                 conf = confFinRModelNew)

# Add fitted values to the data frame
dfFin$fitted_PL <- fitted(FinLModelNew)
dfFin$fitted_PW <- fitted(FinWModelNew)
dfFin$fitted_PR <- fitted(FinRModelNew)

# Plots of fitted values
ggplot(dfFin, aes(x = DOC, y = fitted_PL, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) # fitted, pec fin lengths

ggplot(dfFin, aes(x = DOC, y = fitted_PW, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) # fitted, pec fin widths

ggplot(dfFin, aes(x = DOC, y = fitted_PR, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) # fitted, pec fin ratio (size standardized)

# Plots of non-fitted values
## raw
ggplot(dfFin, aes(x = DOC, y = pecFinLength, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) # raw, pec fin lengths

ggplot(dfFin, aes(x = DOC, y = pecFinBaseWidth, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) # raw, pec fin widths

## size-standardized
ggplot(dfFin, aes(x = DOC, y = finLengthSS, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) # raw, size-standardized pec fin lengths

ggplot(dfFin, aes(x = DOC, y = finBaseSS, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) # raw, size-standardized pec fin widths

ggplot(dfFin, aes(x = DOC, y = finRatioSS, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) # raw, size-standardized pec fin ratio

# Pec Fin Angles ----------------------------------------------------------
dfangle <- read.csv(here("data", "outputs", "PecFinAnglesFINAL.csv")) %>%
  rename("finangle" = "pecFinInsertionAngle")

# Check structure of dfangle
str(dfangle)

# Set basin to factor
dfangle$basin <- as.factor(dfangle$basin)

# With basin
FinAModelNew <- lmer(log(finangle) ~ 1 + log(DOC) + basin + 
                       basin:log(DOC) + (1|lakeID), data = dfangle)
summary(FinAModelNew)
plot(fitted(FinAModelNew) ~ dfangle$DOC, pch = as.numeric(dfangle$basin))

# Calculate profile CIs for parameter estimates
confFinAModelNew <- confint(FinAModelNew, level = 0.95, method = 'profile')

# Assemble parameter estimates, CIs, R2
tableFinAModelNew <- makeSummary(model = FinAModelNew, 
                                 conf = confFinAModelNew)

# Get fitted values and add to data frame. 
dfangle$fitted <- fitted(FinAModelNew)

# Plots:
ggplot(dfangle, aes(x = DOC, y = fitted, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) 

# Raw Plot
ggplot(dfangle, aes(x = DOC, y = finangle, label= lakeID)) + 
  geom_point(aes(colour = lakeID)) 

# Combined univariate model summary tables --------------------------------
tableEyeModelNew <- tableEyeModelNew %>%
  mutate(type = "eyeWidth")
tableRakerLModelNew <- tableRakerLModelNew %>%
  mutate(type = "rakerLength")
tableRakerSModelNew <- tableRakerSModelNew %>%
  mutate(type = "rakerSpace")
tableRakerCModelNew <- tableRakerCModelNew %>%
  mutate(type = "rakerCount")
tableFinLModelNew <- tableFinLModelNew %>%
  mutate(type = "pecFinLength")
tableFinWModelNew <- tableFinWModelNew %>%
  mutate(type = "pecFinWidth")
tableFinRModelNew <- tableFinRModelNew %>%
  mutate(type = "pecFinRatio")
tableFinAModelNew <- tableFinAModelNew %>%
  mutate(type = "pecFinAngle")

# combine into one table
fullTable <- bind_rows(tableEyeModelNew, tableRakerLModelNew, tableRakerSModelNew, tableRakerCModelNew, tableFinLModelNew, tableFinWModelNew, tableFinRModelNew, tableFinAModelNew) %>%
  tidyr::pivot_wider(id_cols = "parameter",
                     names_from = "type",
                     values_from = "estimate")
write.csv(fullTable, here("data", "outputs", "univariateModelSummary.csv"), row.names = F)

# Figures -----------------------------------------------------------------
# The colors and shapes used in these plots are defined in defs.R

# Figures 1, 2, 3 ---------------------------------------------------------
# No code needed for Fig. 1 and Fig. 2--we already have a map and the landmark diagram we can use.

## Fig 3. PC1 vs PC2 plot on DOC gradient, fish body shape. 
### Saved components of this figure above, in sections "PC1 plots" and "PC2 plots"

# Pec Fin Figure ----------------------------------------------------------------
## Pec fins vs. DOC, with gradient. Panels: A (pec fin length), B (pec fin base width), C (length:width ratio), D (insertion anlge)
# Prepare the data, drawing from dfFin and dfangle
df1 <- dfFin %>%
  select(lakeID, DOC, basin, fitted_PL, fitted_PW, fitted_PR) %>%
  distinct()

df2 <- dfangle %>%
  select(lakeID, fitted) %>%
  distinct() %>%
  rename("fitted_A" = fitted)

df <- left_join(df1, df2, by = "lakeID") %>%
  mutate(lakeID = stringr::str_replace(lakeID, "_", " "),
         # Rename basins so they'll show up properly in the legend
         basin = case_when(basin == "4" ~ "Great Lakes basin",
                           basin == "7" ~ "Mississippi basin"),
         # Arrange lake factor levels from high to low DOC
         lakeID = forcats::fct_reorder(lakeID, DOC, .desc = T))

# A (pec fin length)
panelA <- df %>%
  ggplot(aes(x = DOC, y = exp(fitted_PL)))+
  geom_point(aes(fill = lakeID, shape = basin), size = 5, 
             col = "black", stroke = 1)+
  theme_classic()+
  scale_fill_manual(name = "**Lakes**",
                    values = lakeColorsHighLow,
                    guide = guide_legend(override.aes = 
                                           list(shape = lakeShapesHighLow))) +
  scale_shape_manual(name = "",
                     values = c(21, 22))+
  labs(y = "Pectoral fin length (mm)")+
  theme(legend.title = element_markdown(),
        text = element_text(family = "Helvetica", size = 16),
        axis.title.x = element_blank())

# B (pec fin base width)
panelB <- df %>%
  ggplot(aes(x = DOC, y = exp(fitted_PW)))+
  geom_point(aes(fill = lakeID, shape = basin), size = 5, 
             col = "black", stroke = 1)+
  theme_classic()+
  scale_fill_manual(name = "**Lakes**",
                    values = lakeColorsHighLow,
                    guide = guide_legend(override.aes = 
                                           list(shape = lakeShapesHighLow))) +
  scale_shape_manual(name = "",
                     values = c(21, 22))+
  labs(y = "Pectoral fin base width (mm)")+
  theme(legend.title = element_markdown(),
        text = element_text(family = "Helvetica", size = 16),
        axis.title.x = element_blank(),
        legend.position = "none")

# C (length:width ratio)
panelC <- df %>%
  ggplot(aes(x = DOC, y = exp(fitted_PR)))+
  geom_point(aes(fill = lakeID, shape = basin), size = 5, 
             col = "black", stroke = 1)+
  theme_classic()+
  scale_fill_manual(name = "**Lakes**",
                    values = lakeColorsHighLow,
                    guide = guide_legend(override.aes = 
                                           list(shape = lakeShapesHighLow))) +
  scale_shape_manual(name = "",
                     values = c(21, 22))+
  labs(y = "Pectoral fin length:width ratio")+
  theme(legend.title = element_markdown(),
        text = element_text(family = "Helvetica", size = 16),
        axis.title.x = element_blank(),
        legend.position = "none")

# D (insertion angle)
panelD <- df %>%
  ggplot(aes(x = DOC, y = exp(fitted_A)))+
  geom_point(aes(fill = lakeID, shape = basin), size = 5, 
             col = "black", stroke = 1)+
  theme_classic()+
  scale_fill_manual(name = "**Lakes**",
                    values = lakeColorsHighLow,
                    guide = guide_legend(override.aes = 
                                           list(shape = lakeShapesHighLow))) +
  scale_shape_manual(name = "",
                     values = c(21, 22))+
  labs(x = "Dissolved Organic Carbon (mg/L)",
       y = "Pectoral fin insertion angle (degrees)")+
  theme(legend.title = element_markdown(),
        text = element_text(family = "Helvetica", size = 16),
        axis.title.x = element_blank(),
        legend.position = "none")

## extract the legend from panel A
legend <- cowplot::get_legend(panelA)

## make the three stacked plots
allPanels <- cowplot::plot_grid(panelA + 
                                  theme(legend.position = "none"),
                                panelB, panelC, panelD, ncol = 1, align = "v")

## Save the main plot panels
pdf(here("figures", "pecFins", "pecFinsPlotPanels.pdf"), height = 15, width = 5)
allPanels
dev.off()

## Save the legend to the main figures folder
legendPlot <- cowplot::plot_grid(legend)
pdf(here("figures", "legend.pdf"), height = 10, width = 4)
legendPlot
dev.off()

# Make the gradient legend alone
gradientLegend <- cowplot::plot_grid(cowplot::get_legend(df %>%
                                                           ggplot(aes(x = DOC, y = fitted_PL, col = DOC))+
                                                           geom_point()+
                                                           scale_color_gradientn(colors = lakeColorsLowHigh, name = "DOC (mg/L)")))

## Save the gradient legend alone in the main figures folder
pdf(here("figures", "gradientLegend.pdf"), height = 10, width = 4)
gradientLegend
dev.off()

# Gill Raker Figure ----------------------------------------------------------------
## Gill rakers vs. DOC, with gradient. Panels: A (average gill raker length), B (average gill raker space), C (Gill raker total number)
df <- dfraker %>%
  mutate(lakeID = stringr::str_replace(lakeID, "_", " "),
         # Rename basins so they'll show up properly in the legend
         basin = case_when(basin == "4" ~ "Great Lakes basin",
                           basin == "7" ~ "Mississippi basin"),
         # Arrange lake factor levels from high to low DOC
         lakeID = forcats::fct_reorder(lakeID, lakeDOC, .desc = T)) %>%
  select(lakeID, basin, fitted_L, fitted_S, fitted_C, lakeDOC)

## In each plot, we'll exponentiate the fitted values to convert them back to mm units

## Panel A
panelA <- df %>%
  ggplot(aes(x = lakeDOC, y = exp(fitted_L)))+
  geom_point(aes(fill = lakeID, shape = basin), size = 5, 
             col = "black", stroke = 1)+
  theme_classic()+
  scale_fill_manual(name = "**Lakes**",
                    values = lakeColorsHighLow,
                    guide = guide_legend(override.aes = 
                                           list(shape = lakeShapesHighLow))) +
  scale_shape_manual(name = "",
                     values = c(21, 22))+
  labs(y = "Average Gill Raker Length (mm)")+
  theme(legend.title = element_markdown(),
        text = element_text(family = "Helvetica", size = 16),
        axis.title.x = element_blank(),
        legend.position = "none")

## make panel B
panelB <- df %>%
  ggplot(aes(x = lakeDOC, y = exp(fitted_S)))+
  geom_point(aes(fill = lakeID, shape = basin), size = 5, 
             col = "black", stroke = 1)+
  theme_classic()+
  scale_fill_manual(name = "**Lakes**",
                    values = lakeColorsHighLow,
                    guide = guide_legend(override.aes = 
                                           list(shape = lakeShapesHighLow))) +
  scale_shape_manual(name = "",
                     values = c(21, 22))+
  labs(y = "Average Gill Raker Space Width (mm)")+
  theme(legend.title = element_markdown(),
        text = element_text(family = "Helvetica", size = 16),
        axis.title.x = element_blank(),
        legend.position = "none")

## Make panel C, with x-axis label
### not exponentiating the fitted values on this one because they weren't log-transformed.
panelC <- df %>%
  ggplot(aes(x = lakeDOC, y = round(fitted_C, 2)))+
  geom_point(aes(fill = lakeID, shape = basin), size = 5, 
             col = "black", stroke = 1)+
  theme_classic()+
  scale_fill_manual(name = "**Lakes**",
                    values = lakeColorsHighLow,
                    guide = guide_legend(override.aes = 
                                           list(shape = lakeShapesHighLow))) +
  scale_shape_manual(name = "",
                     values = c(21, 22))+
  labs(x = "Dissolved Organic Carbon (mg/L)",
       y = "Gill Raker number")+
  theme(legend.title = element_markdown(),
        text = element_text(family = "Helvetica", size = 16),
        legend.position = "none")

## make the three stacked plots
allPanels <- cowplot::plot_grid(panelA, panelB, panelC, ncol = 1, align = "v")

## Save the main plot panels
pdf(here("figures", "gillRakers", "gillRakersPlotPanels.pdf"), height = 11.25, width = 5)
allPanels
dev.off()

# Eye Width figure ----------------------------------------------------------------
# Eye width across DOC gradient. One panel.
## exponentiating these values to undo the log-transformation
df <- dfeye %>%
  mutate(lakeID = stringr::str_replace(lakeID, "_", " "),
         # Rename basins so they'll show up properly in the legend
         basin = case_when(basin == "4" ~ "Great Lakes basin",
                           basin == "7" ~ "Mississippi basin"),
         # Arrange lake factor levels from high to low DOC
         lakeID = forcats::fct_reorder(lakeID, DOC, .desc = T)) %>%
  select(lakeID, DOC, basin, fitted) %>%
  distinct()

eyePlot <- df %>%
  ggplot(aes(x = DOC, y = exp(fitted)))+
  geom_point(aes(fill = lakeID, shape = basin), size = 5, 
             col = "black", stroke = 1)+
  theme_classic()+
  scale_fill_manual(name = "**Lakes**",
                    values = lakeColorsHighLow,
                    guide = guide_legend(override.aes = 
                                           list(shape = lakeShapesHighLow))) +
  scale_shape_manual(name = "",
                     values = c(21, 22))+
  labs(x = "Dissolved Organic Carbon (mg/L)",
       y = "Eye Width (mm)")+
  theme(legend.title = element_markdown(),
        text = element_text(family = "Helvetica", size = 16),
        legend.position = "none")

## Save the main plot panel
pdf(here("figures", "eyeWidths", "eyeWidthsPlotPanel.pdf"), 
    height = 3.75, width = 5)
eyePlot
dev.off()
