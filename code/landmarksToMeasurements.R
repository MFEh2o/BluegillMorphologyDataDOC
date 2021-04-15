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
db <- "MFEdb_20210402.db" # name of db

# Load the FISH_MORPHOMETRICS data ----------------------------------------
fm <- dbTable("fish_morphometrics")

# Define a function to compute Euclidean distance -------------------------
# There is probably already a function for this, but I didn't feel like loading extra packages and it was simple enough to write!
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
# Distance between landmark 13 and landmark 14 (defined in Chelsea's thesis)
fm %>%
  group_by(imageFile) %>%
  filter("pecFinBaseWidth" %in% parameter) %>%
  pull(parameter) %>% 
  unique()
# XXX question here: even though landmarks 13 and 14 are marked on the diagram of the pec fins alone, we don't have data for landmarks 13 and 14 for the pec-fin-only photos. Because of that, where did the pecFinBaseWidth measurements come from? Did they come from the full-body photos, or were they taken with the ruler in imageJ, or were landmarks 13 and 14 measured on the pec-fin-only photos but not recorded?

# pecFinInsertionAngle ----------------------------------------------------
# Angle between landmark 13 and landmark 14
pfa <- fm %>%
  filter(parameter %in% c("pecFinInsertionAngle", "X13", "Y13", "X14", "Y14")) %>%
  select(imageFile, replicate, parameter, parameterValue) %>%
  pivot_wider(id_cols = c("imageFile", "replicate"), names_from = parameter, values_from = parameterValue) %>%
  filter(!is.na(X13), !is.na(X14))

pfa$pecFinInsertionAngleRecomputed <- NA # initialize column
## Attempt to use angle.calc function from Morpho package to compute angles
for(i in 1:nrow(pfa)){
  pfa$pecFinInsertionAngleRecomputed[i] <- angle.calc(x = c(pfa$X14[i], pfa$Y14[i]),
                                                  y = c(pfa$X13[i], pfa$Y13[i]))*(180/pi) # divide by 180/pi to convert from radians to degrees

} # XXX THIS ISN'T WORKING--angles don't make sense. Where could I find the original angle calculation?


# pecFinLength ------------------------------------------------------------
# Distance between landmark 14 and landmark 21
# XXX same question here as for the pecFinBaseWidth.
fm %>%
  group_by(imageFile) %>%
  filter("pecFinLength" %in% parameter) %>%
  pull(parameter) %>% 
  unique()


# pecFinRatioSizeStd: "width of the insertion point of the fin divided by its total length." (size-standardized)
# raker lengths: measured with the ruler in imageJ or tpsDig
# raker spaces: measured with the ruler in imageJ or tpsDig


