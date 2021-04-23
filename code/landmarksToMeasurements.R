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

pfa <- pfa %>%
  mutate(pecFinInsertionAngleRecomputed = 
           asin( # arcsin(opposite/hypotenuse) = angle (see diagram in Resources/ folder)
             (Y14-Y13)/ # Vertical distance between 14 and 14 (opposite)
               eudist(X13, X14, Y13, Y14) # eudist function defined above to compute linear distance between two points. Arguments in the following order: X1, X2, Y1, Y2. (hypotenuse)
             )*
           (180/pi) # multiply by 180/pi to convert from radians to degrees
         )

head(pfa %>% select(contains("pecFinInsertionAngle")), 10) %>% mutate(sum = pecFinInsertionAngle + pecFinInsertionAngleRecomputed) # wow, this does not match up at all with pecFinInsertionAngle!


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


