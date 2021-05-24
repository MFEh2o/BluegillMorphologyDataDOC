# Main analysis script
# Created by Chelsea Bishop. Edits/reorganization by Kaija Gahm
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
## landmark data
myData_tps <- readland.tps(here("data", "inputs", "FULL_2018_TPS_FILE_UPDATED_09-25-19.TPS"), specID = "imageID")

## manually-defined link data
mylinks <- read.table(here("data", "inputs", "Full_body_links.txt"))

## fish metadata to correspond to the landmarks
identifiers <- read.table(here("data", "outputs", "Identifiers_Update_2020.txt"), 
                          sep = ",", header = TRUE)
# Landmark repeatability analysis -----------------------------------------
# Following the procedure outlined on page 258 of Zelditch 2012, and in this blog post: https://www.r-bloggers.com/2015/04/tips-tricks-8-examining-replicate-error/
# 1. Read the coordinate data into R (and format it)
tps_rep <- readland.tps(here("data", "inputs", "Total Replicates Coords.TPS"),
                        specID = "imageID")

### Make data frame
gdf.rep <- geomorph.data.frame(shape = tps_rep)

# 2. Use gpagen() to perform a Procrustes Superimposition
superimp <- gpagen(tps_rep, ProcD = TRUE, Proj = TRUE)  
plot(superimp) # all of these look reasonable--no landmarks way far away from where they should be.
corePCA <- plotTangentSpace(superimp$coords)

### Classifiers
classifiers <- read.table(here("data", "outputs", "Replicates Identifiers.txt"), 
                          sep = ",", header = TRUE)
gdf.rep$replic <- as.factor(classifiers$replicateNo)
gdf.rep <- geomorph.data.frame(shape = tps_rep, 
                               replic = classifiers$replicateNo, 
                               fishID = classifiers$imageID,
                               cSize = superimp$Csize)

# 3. Perform a Procrustes ANOVA
### ANOVA
fishrep.anova <- procD.lm(shape ~ fishID + fishID:replic, data = gdf.rep) 
summary(fishrep.anova)
# XXX MADLEN SAYS ONE OF THESE HAS TO BE SIGNIF AND THE OTHER NOT. BUT I SEE BOTH SIGNIF?
# XXX GRANT DISAGREES WITH MADLEN.

###             Df     SS      MS     Rsq        F       Z Pr(>SS)   
#fishID         59 602102 10205.1 0.97561 187.1758 10.3693   0.001 **
#fishID:replic  60   8509   141.8 0.01379   2.6012  5.3186   0.001 **
#Residuals     120   6543    54.5 0.01060                            
#Total         239 617154        

# To calculate the repeatability of our digitizing ability, we subtract the MS of the replicate term from the individual term and divide by 4 (because we measured each fish 4 times):
#  (MS(fishID) – MS(fishID:replic))/4
# To make these calculations easier, I'm going to extract the MS terms of the ANOVA output as their own variables with descriptive names: 
ms_inter <- fishrep.anova$aov.table[1,3] # MS term for the among-fish variability
ms_intra <- fishrep.anova$aov.table[2,3] # MS term for the within-individual, among-replicates variability
part1 <- (ms_inter - ms_intra)/4

# Then we calculate the ratio of this value (part1) to the total MS:
# ((MS(fishID) – MS(fishID:replic))/4 ) / (MS(fishID) + MS(fishID:replic)) 

repeatability <- part1/(ms_intra + part1)
# repeatability is 0.947. # XXX DIFFERS FROM THE MANUSCRIPT
# Note that the blog post referenced above would have us doing repeatability <- part/(ms_intra + ms_inter), but that gives a value around 0.24, which just *can't* possibly be right. I'm going with Zelditch instead.
# Zelditch says "measurement error is often quantified as repeatability (R) using a ratio of two variance components, that for among-individual to the sum of the among-individual and measurement error components.
# The terminology there is really confusing, and Zelditch provides methods of calculating each of these components that I can't quite follow. But by combining this description with the method outlined in the blog post above for extracting each MS from the ANOVA, we do end up with R = part1/(ms_intra + part1)

# Prepare the data for morphometric analyses -------------------------------
## set imageID's as names for the individual fish (third dimension)
dimnames(myData_tps)[[3]] <- identifiers$imageID

## number the landmarks (first dimension)
dimnames(myData_tps)[[1]] <- as.character(1:19)

# Morphometric analysis ------------------------------------------------------
## perform Procrustes analysis
GPA.fish <- gpagen(myData_tps, ProcD = TRUE, Proj = TRUE) 

## convert results of Procrustes analysis into a usable data frame
gdf.fish <- geomorph.data.frame(shape = GPA.fish$coords,
                                Lake = identifiers$lakeID,
                                DOC = identifiers$lakeDOC,
                                captureMethod = identifiers$captureMethod, 
                                cSize = GPA.fish$Csize, 
                                basin = as.factor(identifiers$Basin))

## PCA scoreplot for overall fish shapes, colored by DOC levels (proxy for lakes)
corePCA <- plotTangentSpace(gdf.fish$shape, 
                            groups = as.factor(gdf.fish$DOC), 
                            legend = TRUE)

## save PC scores to their own object
PCscores <- corePCA$pc.scores[,1:2]

## Split PC scores into individual objects by lake
Bay <- PCscores[which(gdf.fish$Lake == "BA"),]
Birch <- PCscores[which(gdf.fish$Lake == "BH"),]
Crampton <- PCscores[which(gdf.fish$Lake == "CR"),]
Found <- PCscores[which(gdf.fish$Lake == "FD"),]
Hummingbird <- PCscores[which(gdf.fish$Lake == "HB"),]
Little_Crooked <- PCscores[which(gdf.fish$Lake == "LC"),]
Lost <- PCscores[which(gdf.fish$Lake == "LT"),]
McCullough <- PCscores[which(gdf.fish$Lake == "MC"),]
Muskellunge <- PCscores[which(gdf.fish$Lake == "MS"),]
Oxbow <- PCscores[which(gdf.fish$Lake == "OB"),]
Papoose <- PCscores[which(gdf.fish$Lake == "PS"),]
Red_Bass <- PCscores[which(gdf.fish$Lake == "RS"),]
Squaw <- PCscores[which(gdf.fish$Lake == "SQ"),]
Towanda <- PCscores[which(gdf.fish$Lake == "TO"),]


# Save plots of the PC1 extremes ------------------------------------------
# Open a pdf file
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

# Open a pdf file
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
# Open a pdf file
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

# Open a pdf file
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
gpa_array <- gpagen(myData_tps, ProcD = TRUE, Proj = TRUE)$coords # takes just the coords

# Convert array to matrix for PCA
gpa_mat <- t(apply(gpa_array, 3, function(y) matrix(t(y), 1)))

# Perform non-phylogenetic PCA
resEig <- eigen(cov(gpa_mat))

# Get PC scores, changed sign of rotation matrix (resEig$vectors)
scores <- gpa_mat %*% -(resEig$vectors)
#### shapes switched, need to reverse signs to fix

# Get percent variance explained along each axis
per_var <- (resEig$values / sum(resEig$values))*100

# Set PCs to plot
pcs <- 1:2

# (all individuals) Backtransform morphospace all individuals -------------------------------
# NOTE: Because I wasn't sure what Chelsea would want to use for her photoshopping, I've created three versions of this plot, all on the axes defined from a PCA on individual points (not lake means). One plot shows the lake means unlabeled, another shows the labeled lake means, and the last shows all individual points (in case that's desirable/useful for something).

# Create data frame of means by lake, for use later in plotting the legend
means <- data.frame(mn1 = c(mean(Hummingbird[,1]),
                            mean(Squaw[,1]),
                            mean(Red_Bass[,1]),
                            mean(McCullough[,1]),
                            mean(Oxbow[,1]),
                            mean(Birch[,1]),
                            mean(Bay[,1]),
                            mean(Muskellunge[,1]),
                            mean(Papoose[,1]),
                            mean(Found[,1]),
                            mean(Towanda[,1]),
                            mean(Crampton[,1]),
                            mean(Little_Crooked[,1]),
                            mean(Lost[,1])),
                    mn2 = c(mean(Hummingbird[,2]),
                            mean(Squaw[,2]),
                            mean(Red_Bass[,2]),
                            mean(McCullough[,2]),
                            mean(Oxbow[,2]),
                            mean(Birch[,2]),
                            mean(Bay[,2]),
                            mean(Muskellunge[,2]),
                            mean(Papoose[,2]),
                            mean(Found[,2]),
                            mean(Towanda[,2]),
                            mean(Crampton[,2]),
                            mean(Little_Crooked[,2]),
                            mean(Lost[,2])))

# 1. PLOT WITH MEAN POINTS ONLY
# Open a pdf file
pdf(here("figures", "fishShapes_pc1_pc2", "mainPlot_meansOnly.pdf"), width = 9.5, height = 7) 
# Create plot box with axes and axis labels
par(mar = c(6, 6, 6, 6)) # set margins around the plot so we can read the axis labels
plot(scores[, pcs], type = "n", # type = "n" plots the axes without plotting the scores, so we can add the scores later
     main = "Backtransform morphospace", # title
     xlab = paste0("PC", pcs[1], " (", round(per_var[pcs[1]]), "%)"),
     ylab = paste0("PC", pcs[2], " (", round(per_var[pcs[2]]), "%)"),
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5) # the cex. commands scale up the text to make it larger and more easily readable

# Plot backtransform shapes, changed sign of rotation matrix (resEig$vectors) 
# There's no documentation for the btShapes function (i.e. ?btShapes or help(btShapes) won't get you anything), but see this tutorial (https://aaronolsen.github.io/tutorials/morphometrics/backtransform.html) for an explanation of the input parameters to btShapes. 
StereoMorph::btShapes(scores = scores, 
                      vectors = -(resEig$vectors), 
                      fcn = plot_fish_lateral, # defined in plotFishLateral.R
         pcs = pcs, # the two pc axes for the plot
         n = c(4,4), # number of shapes to draw along the x and y axes, respectively
         m = dim(lm_array)[2], 
         row.names= dimnames(lm_array)[[1]], # order in which to connect the landmarks
         pc.margin = c(0.06, 0.05), size = 0.038, col = gray(0.7))

# add points for each lake in a different color
points(mean(Hummingbird[,1]), mean(Hummingbird[,2]), col=lakeColorsHighLow[1], pch=19,cex=3)
points(mean(Squaw[,1]), mean(Squaw[,2]), col=lakeColorsHighLow[2], pch=19,cex=3)
points(mean(Red_Bass[,1]), mean(Red_Bass[,2]), col=lakeColorsHighLow[3], pch=19,cex=3)
points(mean(McCullough[,1]), mean(McCullough[,2]), col=lakeColorsHighLow[4], pch=19,cex=3)
points(mean(Oxbow[,1]), mean(Oxbow[,2]), col=lakeColorsHighLow[5], pch=19,cex=3)
points(mean(Birch[,1]), mean(Birch[,2]), col=lakeColorsHighLow[6], pch=19,cex=3)
points(mean(Bay[,1]), mean(Bay[,2]), col=lakeColorsHighLow[7], pch=19,cex=3)
points(mean(Muskellunge[,1]), mean(Muskellunge[,2]), col=lakeColorsHighLow[8], pch=19,cex=3)
points(mean(Papoose[,1]), mean(Papoose[,2]), col=lakeColorsHighLow[9], pch=19,cex=3)
points(mean(Found[,1]), mean(Found[,2]), col=lakeColorsHighLow[10], pch=19,cex=3)
points(mean(Towanda[,1]), mean(Towanda[,2]), col=lakeColorsHighLow[11], pch=19,cex=3)
points(mean(Crampton[,1]), mean(Crampton[,2]), col=lakeColorsHighLow[12], pch=19,cex=3)
points(mean(Little_Crooked[,1]), mean(Little_Crooked[,2]), col=lakeColorsHighLow[13], pch=19,cex=3)
points(mean(Lost[,1]), mean(Lost[,2]), col=lakeColorsHighLow[14], pch=19,cex=3)

legend("topright", legend = lakesHighLow, 
       pch = 19, col = lakeColorsHighLow)
dev.off()

# 2. PLOT WITH LABELED MEAN POINTS
# Open a pdf file
pdf(here("figures", "fishShapes_pc1_pc2", "mainPlot_labeledMeans.pdf"), width = 9.5, height = 7) 
# Create plot box with axes and axis labels
par(mar = c(6, 6, 6, 6)) 
plot(scores[, pcs], type = "n", main = "Backtransform morphospace",
     xlab = paste0("PC", pcs[1], " (", round(per_var[pcs[1]]), "%)"),
     ylab = paste0("PC", pcs[2], " (", round(per_var[pcs[2]]), "%)"),
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)

# Plot backtransform shapes, changed sign of rotation matrix (resEig$vectors) 
StereoMorph::btShapes(scores = scores, vectors = -(resEig$vectors), fcn = plot_fish_lateral, 
                      pcs = pcs, n = c(4,4), m = dim(lm_array)[2], 
                      row.names= dimnames(lm_array)[[1]], 
                      pc.margin = c(0.06,0.05), size = 0.038, col = gray(0.7))

# add points for each lake in a different color
points(mean(Hummingbird[,1]), mean(Hummingbird[,2]), col=lakeColorsHighLow[1], pch=19,cex=3)
points(mean(Squaw[,1]), mean(Squaw[,2]), col=lakeColorsHighLow[2], pch=19,cex=3)
points(mean(Red_Bass[,1]), mean(Red_Bass[,2]), col=lakeColorsHighLow[3], pch=19,cex=3)
points(mean(McCullough[,1]), mean(McCullough[,2]), col=lakeColorsHighLow[4], pch=19,cex=3)
points(mean(Oxbow[,1]), mean(Oxbow[,2]), col=lakeColorsHighLow[5], pch=19,cex=3)
points(mean(Birch[,1]), mean(Birch[,2]), col=lakeColorsHighLow[6], pch=19,cex=3)
points(mean(Bay[,1]), mean(Bay[,2]), col=lakeColorsHighLow[7], pch=19,cex=3)
points(mean(Muskellunge[,1]), mean(Muskellunge[,2]), col=lakeColorsHighLow[8], pch=19,cex=3)
points(mean(Papoose[,1]), mean(Papoose[,2]), col=lakeColorsHighLow[9], pch=19,cex=3)
points(mean(Found[,1]), mean(Found[,2]), col=lakeColorsHighLow[10], pch=19,cex=3)
points(mean(Towanda[,1]), mean(Towanda[,2]), col=lakeColorsHighLow[11], pch=19,cex=3)
points(mean(Crampton[,1]), mean(Crampton[,2]), col=lakeColorsHighLow[12], pch=19,cex=3)
points(mean(Little_Crooked[,1]), mean(Little_Crooked[,2]), col=lakeColorsHighLow[13], pch=19,cex=3)
points(mean(Lost[,1]), mean(Lost[,2]), col=lakeColorsHighLow[14], pch=19,cex=3)

legend("topright", legend = lakesHighLow, 
       pch = 19, col = lakeColorsHighLow)
text(means[,1], means[,2], labels = lakesHighLow, cex = 1.5)
dev.off()

# 3. PLOT WITH all points
# Open a pdf file
pdf(here("figures", "fishShapes_pc1_pc2", "mainPlot_allIndividuals.pdf"), width = 9.5, height = 7) 
# Create plot box with axes and axis labels
par(mar = c(6, 6, 6, 6)) 
plot(scores[, pcs], type = "n", main = "Backtransform morphospace",
     xlab = paste0("PC", pcs[1], " (", round(per_var[pcs[1]]), "%)"),
     ylab = paste0("PC", pcs[2], " (", round(per_var[pcs[2]]), "%)"),
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)

# Plot backtransform shapes, changed sign of rotation matrix (resEig$vectors) 
StereoMorph::btShapes(scores = scores, vectors = -(resEig$vectors), fcn = plot_fish_lateral, 
                      pcs = pcs, n = c(4,4), m = dim(lm_array)[2], 
                      row.names= dimnames(lm_array)[[1]], 
                      pc.margin = c(0.06,0.05), size = 0.038, col = gray(0.7))

# add points for each lake in a different color
points(Hummingbird[,1], Hummingbird[,2], col=lakeColorsHighLow[1], pch=19,cex=2,alpha = 0.7)
points(Squaw[,1], Squaw[,2], col=lakeColorsHighLow[2], pch=19,cex=2,alpha = 0.7)
points(Red_Bass[,1], Red_Bass[,2], col=lakeColorsHighLow[3], pch=19,cex=2,alpha = 0.7)
points(McCullough[,1], McCullough[,2], col=lakeColorsHighLow[4], pch=19,cex=2,alpha = 0.7)
points(Oxbow[,1], Oxbow[,2], col=lakeColorsHighLow[5], pch=19,cex=2,alpha = 0.7)
points(Birch[,1], Birch[,2], col=lakeColorsHighLow[6], pch=19,cex=2,alpha = 0.7)
points(Bay[,1], Bay[,2], col=lakeColorsHighLow[7], pch=19,cex=2,alpha = 0.7)
points(Muskellunge[,1], Muskellunge[,2], col=lakeColorsHighLow[8], pch=19,cex=2,alpha = 0.7)
points(Papoose[,1], Papoose[,2], col=lakeColorsHighLow[9], pch=19,cex=2,alpha = 0.7)
points(Found[,1], Found[,2], col=lakeColorsHighLow[10], pch=19,cex=2,alpha = 0.7)
points(Towanda[,1], Towanda[,2], col=lakeColorsHighLow[11], pch=19,cex=2,alpha = 0.7)
points(Crampton[,1], Crampton[,2], col=lakeColorsHighLow[12], pch=19,cex=2,alpha = 0.7)
points(Little_Crooked[,1], Little_Crooked[,2], col=lakeColorsHighLow[13], pch=19,cex=2,alpha = 0.7)
points(Lost[,1], Lost[,2], col=lakeColorsHighLow[14], pch=19,cex=2,alpha = 0.7)

legend("topright", legend = lakesHighLow, 
       pch = 19, col = lakeColorsHighLow)
dev.off()

# (per-lake means) Backtransform morphospace figure with mean shapes per lake ----------
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

# Plot
# Create plot box with axes and axis labels
par(mar = c(6, 6, 6, 6)) 
plot(scores[, pcs], type = "n", main = "Backtransform morphospace",
     xlab = paste0("PC", pcs[1], " (", round(per_var[pcs[1]]), "%)"),
     ylab = paste0("PC", pcs[2], " (", round(per_var[pcs[2]]), "%)"), 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)


# Plot backtransform shapes, changed sign of rotation matrix (resEig$vectors) 
btShapes(scores = scores, vectors = -(resEig$vectors), 
         fcn = plot_fish_lateral, 
         pcs = pcs, n = c(4,4), m = dim(lm_array)[2], 
         row.names = dimnames(lm_array)[[1]], 
         pc.margin = c(0.06,0.06), 
         size = 0.018, 
         col = gray(0.7))
points(scores [,1], scores[,2], col = "#41C6EC", pch = 19, cex = 2)
text(scores [,1], scores [,2], labels = levels(gpa_df$Lake), cex = 1.5)
# For ref-to-Target Figures that go along the axes in figure 3: see pc1Plot and pc2Plot, above.

# Stats for shape ---------------------------------------------------------

#Model structure similar to univariate analyses

#Check that basin is a factor
str(gdf.fish)

# Fit model
# shapeModel <- procD.lm(shape ~ log(cSize) + log(DOC) + basin + log(DOC):basin + Lake, data = gdf.fish, SS.type = "II")

# Because this is not a hierarchical model, the log(DOC):basin term is redundant with the Lake term. Drop the former.
shapeModel <- procD.lm(shape ~ log(cSize) + log(DOC) + basin + Lake, 
                       data = gdf.fish, SS.type = "II", iter = 10000)

#Diagnostic plot
plot(shapeModel)

#ANOVA
anova(shapeModel,error = c("Residuals", "Residuals", "Residuals", "Lake"))

# Univariate Data Check and Run -------------------------------------------
# Eye Widths --------------------------------------------------------------
dfeye <-read.csv(here("data", "outputs", "eyewidthsFINAL.csv")) 

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
  geom_point(aes(colour = lakeID))
ggplot(dfeye,aes(x = DOC, y = eyewidth, label = lakeID)) + 
  geom_point(aes(colour = lakeID))
ggplot(dfeye,aes(x = DOC, y = eyewidth.ss, label = lakeID)) + 
  geom_point(aes(colour = lakeID))

# Gill Rakers -------------------------------------------------------------
dfraker <- read.csv(here("data", "outputs", "Gill_Rakers_2018_Final.csv"))

# Check structure of dfraker
str(dfraker)

# Set basin to factor
dfraker$basin <- as.factor(dfraker$basin)

# avgL2_ss is size standardizations for average length for rakers 4-7

## lengths
RakerLModelNew <- lmer(log(avgL2_ss) ~ 1 + log(lakeDOC) + basin + 
                         basin:log(lakeDOC) + (1|lakeID), data = dfraker)
summary(RakerLModelNew)

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
nrow(dfFin) # still 218

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
  geom_point(aes(colour = lakeID)) 

ggplot(dfFin, aes(x = DOC, y = fitted_PW, label = lakeID)) + 
  geom_point(aes(colour = lakeID))

ggplot(dfFin, aes(x = DOC, y = fitted_PR, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) 

# Plots of non-fitted values
## raw
ggplot(dfFin, aes(x = DOC, y = pecFinLength, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) 

ggplot(dfFin, aes(x = DOC, y = pecFinBaseWidth, label = lakeID)) + 
  geom_point(aes(colour = lakeID))

## size-standardized
ggplot(dfFin, aes(x = DOC, y = finLengthSS, label = lakeID)) + 
  geom_point(aes(colour = lakeID)) 

ggplot(dfFin, aes(x = DOC, y = finBaseSS, label = lakeID)) + 
  geom_point(aes(colour = lakeID))

ggplot(dfFin, aes(x = DOC, y = finRatioSS, label = lakeID)) + 
  geom_point(aes(colour = lakeID))

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