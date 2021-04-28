# Script to double-check that the measurement data in FISH_MORPHOMETRICS can be re-created from the landmarks
# Created by Kaija Gahm on 7 April 2021

# Load packages -----------------------------------------------------------
library(dplyr) # for data wrangling
library(tidyr) # for data wrangling
library(stringr) # for text string manipulation
library(here) # for file paths
library(ggplot2)
source(here("code", "dbUtil.R")) # for the dbTable() function
library(RSQLite) # for db connection
library(Morpho) # for angles

# Connect to the database -------------------------------------------------
dbdir <- here("data") # directory where the db is stored
db <- "MFEdb_20210423.db" # name of db

# Load the FISH_MORPHOMETRICS data ----------------------------------------
fm <- dbTable("fish_morphometrics")

# Define a function to compute Euclidean distance -------------------------
eudist <- function(x1, x2, y1, y2){ # takes four vectors as input
  deltax <- x2-x1
  deltay <- y2-y1
  dist <- sqrt(deltax^2 + deltay^2)
  return(dist)
}

# eyeWidth ----------------------------------------------------------------
# Distance between landmark 2 and landmark 19 (from Chelsea's thesis)
ew <- fm %>%
  filter(parameter %in% c("eyeWidth", "X2", "Y2", "X19", "Y19")) %>%
  select(imageFile, replicate, parameter, parameterValue) %>%
  pivot_wider(id_cols = c("imageFile", "replicate"), names_from = parameter, values_from = parameterValue) %>%
  filter(!is.na(X2), !is.na(X19))

ew <- ew %>%
  mutate(eyeWidthRecomputed = eudist(X2, X19, Y2, Y19), # recompute eyewidth
         diff = eyeWidth - eyeWidthRecomputed) # difference between the recomputed values and the original values

## Plot the 'diff' results to see if there are any concerning differences
ew %>%
  ggplot(aes(x = diff))+
  geom_density() # looks good--tiny differences due to rounding, but we don't need to worry about that.


# pecFinBaseWidth ---------------------------------------------------------
# Distance between landmark 13 and landmark 14, computed on PEC FIN photos
recomp <- fm %>%
  filter(imageType == "pec_fin_only",
         parameter %in% c("pecFinBaseWidth", "X13", "Y13", "X14", "Y14")) %>%
  pivot_wider(id_cols = c("imageFile", "replicate"),
              names_from = "parameter",
              values_from = "parameterValue") %>%
  mutate(pecFinBaseWidthRecomputed = eudist(X13, X14, Y13, Y14),
         diff = pecFinBaseWidthRecomputed - pecFinBaseWidth)

# Visualize the differences
recomp %>%
  ggplot(aes(x = diff))+
  geom_density() # yeah, these are all basically zero--rounding differences. Not a concern. Good!

# pecFinInsertionAngle ----------------------------------------------------
# I'm not going to go through the steps of recomputing this because I just computed the angles now, and I'd just be rehashing my old code. 

# "pecFinAngleDiagram.png" in Resources/ shows how this angle was computed--it's the angle marked in red.
# The formula is (code copied from "gh142_pecFinPhotoLandmarks.R" for the MFE database--version 4.5.4, 2021-04-23:

# angle = 
#   asin( # arcsin(opposite/hypotenuse) = angle (see diagram in Resources/ folder)
#     (Y14-Y13)/ # Vertical distance between 14 and 14 (opposite)
#       eudist(X13, X14, Y13, Y14) # eudist function defined above to compute linear distance between two points. Arguments in the following order: X1, X2, Y1, Y2. (hypotenuse)
#   )*
#   (180/pi) # multiply by 180/pi to convert from radians to degrees

# The function "eudist" is defined above. Euclidean distance between landmarks 13 and 14, in this case.

# pecFinLength ------------------------------------------------------------
# Distance between landmark 14 and landmark 21
# Likewise, I'm not going to bother here, because I recomputed these lengths already and found some discrepancies. I used that to make updates to FISH_MORPHOMETRICS in the gh142_pecFinPhotoLandmarks.R script in the MFE database pipeline. So they've been recently checked.

# pecFinRatioSizeStd ------------------------------------------------------
  


# pecFinRatioSizeStd: "width of the insertion point of the fin divided by its total length." (size-standardized)


# raker lengths: measured with the ruler in imageJ or tpsDig
# raker spaces: measured with the ruler in imageJ or tpsDig


