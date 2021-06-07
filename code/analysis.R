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
# repeatability is 0.947. # XXX DIFFERS FROM THE MANUSCRIPT
# Note that the blog post referenced above would have us doing repeatability <- part/(ms_intra + ms_inter), but that gives a value around 0.24, which just *can't* possibly be right. I'm going with Zelditch instead.
# Zelditch says "measurement error is often quantified as repeatability (R) using a ratio of two variance components, that for among-individual to the sum of the among-individual and measurement error components.
# The terminology there is really confusing, and Zelditch provides methods of calculating each of these components that I can't quite follow. But by combining this description with the method outlined in the blog post above for extracting each MS from the ANOVA, we do end up with R = part1/(ms_intra + part1)

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
# NOTE: Because I wasn't sure what Chelsea would want to use for her photoshopping, I've created three versions of this plot, all on the axes defined from a PCA on individual points (not lake means). One plot shows the lake means unlabeled, another shows the labeled lake means, and the last shows all individual points (in case that's desirable/useful for something).

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
# Wrapper function is defined in defs.R
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

# Plot backtransform shapes, changed sign of rotation matrix (resEig_indiv$vectors) 
# There's no documentation for the btShapes function (i.e. ?btShapes or help(btShapes) won't get you anything), but see this tutorial (https://aaronolsen.github.io/tutorials/morphometrics/backtransform.html) for an explanation of the input parameters to btShapes. 
# Wrapper function is defined in defs.R
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
#### shapes switched need to reverse signs to fix

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
  mag=1, 
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
  # target = min (orange)
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
  # target = max (orange)
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
  mag=1, 
  useRefPts = TRUE)
dev.off()

# XXX STILL NEED TO ADD THE REST OF THE ANALYSIS
