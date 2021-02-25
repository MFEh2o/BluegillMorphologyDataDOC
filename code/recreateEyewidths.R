# Script to recreate eyewidthsFINAL.csv
# Created by Kaija Gahm on 19 February 2021

# Load packages -----------------------------------------------------------
library(dplyr) # for data wrangling
library(tidyr) # for data wrangling
library(stringr) # for text string manipulation
library(here) # for file paths
source(here("code", "dbUtil.R")) # for the dbTable() function
library(RSQLite) # for db connection

# Connect to the database -------------------------------------------------
dbdir <- here("data") # directory where the db is stored
db <- "MFEdb_20201125.db" # name of db

# Load the original eyewidthsFINAL for comparison -------------------------
ew <- read.csv(here("data", "unclassified", "eyewidthsFINAL_ORIGINAL.csv")) %>%
  select(-imageID) # this is a file path to chelsea's computer, so don't need this.

# Grab FISH_MORPHOMETRICS -------------------------------------------------
fm <- dbTable("fish_morphometrics")

# Grab lakeInfo -----------------------------------------------------------
lakeInfo <- read.csv(here("data", "outputs", "lakeInfo_wBins.csv"))

# Initialize recreated df -------------------------------------------------
ewR <- fm %>%
  select(imageFile, parameter, parameterValue) %>%
  rename("fishID" = imageFile) %>%
  filter(parameter == "eyeWidth") %>%
  pivot_wider(id_cols = fishID, names_from = "parameter", values_from = "parameterValue")

# Are all the fish represented?
all(ew$fishID %in% ewR$fishID) # good!

# Limit it to the fish contained in the original file
ewR <- ewR %>%
  filter(fishID %in% ew$fishID)

nrow(ewR) == nrow(ew)
all(ewR$fishID == ew$fishID)
all(ewR$eyeWidth == ew$eyeWidth) # all right, these match!

# lakeID ------------------------------------------------------------------
ewR <- ewR %>%
  mutate(lakeID = factor(str_extract(fishID, "^[A-Z]+(?=\\s)"))) %>%
  mutate(lakeID = forcats::fct_recode(lakeID,
                                      "Bay" = "BA",
                                      "Birch" = "BH",
                                      "Crampton" = "CR",
                                      "Found" = "FD",
                                      "Hummingbird" = "HB",
                                      "Little_Crooked" = "LC",
                                      "Lost" = "LT",
                                      "McCullough" = "MC",
                                      "Muskellunge" = "MK",
                                      "Oxbow" = "OB",
                                      "Papoose" = "PP",
                                      "Red_Bass" = "RB",
                                      "Squaw" = "SQ",
                                      "Towanda" = "TW"))

table(ewR$lakeID, exclude = NULL)
table(ew$lakeID, exclude = NULL) # good, the counts line up and there are no NA's. I'm using lake names to avoid incorrect abbreviations.

# Add DOC and basin ------------------------------------------------------
ewR <- ewR %>%
  left_join(lakeInfo %>%
              select(lakeName, DOC, DOClevel, basin,
                     "area" = surface_Area.ha., 
                     "max_Depth" = max_Depth.m.),
            by = c("lakeID" = "lakeName"))

## Check that the DOC values went through
all(ewR$DOC == ew$DOC) # yay!

# Not adding capture method because it doesn't get used in analyses.

# Calculate fish standard length ------------------------------------------
# Standard length is the distance between landmarks 1 and 7
fsl <- fm %>%
  filter(parameter %in% c("X7", "Y7", "X1", "Y1")) %>%
  select(imageFile, parameter, parameterValue) %>%
  group_by(imageFile, parameter) %>% # because standard length was calculated based on just the first measurement, we group by imageFile and parameter and select the first row.
  slice(1) %>%
  ungroup() %>%
  pivot_wider(names_from = parameter, values_from = parameterValue) %>%
  rename("fishID" = imageFile) %>%
  mutate(fishStdLength = sqrt((X7-X1)^2+(Y7-Y1)^2)) %>%
  select(fishID, fishStdLength) %>%
  filter(fishID %in% ewR$fishID)

## join
ewR <- ewR %>%
  left_join(fsl, by = "fishID")

## Check
ewR$fishStdLength %>% head(10)
ew$fishStdLength %>% head(10)

# ok, these are the same to 2 decimal places. So not *quite* exactly the same. But it looks like this column actually doesn't get used in the modeling, so that's okay.

# Size-standardize the eye widths -----------------------------------------
## Models to account for/remove fish size
## First, a model that includes an interaction effect with `lakeID`:
size.rem <- lm(log10(eyeWidth) ~ log10(fishStdLength)*lakeID, data = ewR)
summary(size.rem) # we note that there is an interaction effect for CR, LT, and MS.

## To get the common within-group slope to be used in size-standardization, we need to remove that interaction effect. So we'll fit a second model:
size.rem2 <- lm(log10(eyeWidth) ~ log10(fishStdLength)+lakeID, data = ewR)
summary(size.rem2) # now we can extract the common within-group slope, which is the `Estimate` parameter for `log10(fishStdLength)`: 0.626646.

commonWithinGroupSlope <- coef(size.rem2)[2] %>% unname() # programmatically grab that coefficient
meanFishSize <- mean(ewR$fishStdLength) # compute the mean `fishStdLength` across all fish in the data set

## Now add a column to ewR for size-standardized eye width
ewR <- ewR %>%
  mutate(eyewidth.ss = eyeWidth*(meanFishSize/fishStdLength)^commonWithinGroupSlope)

## check it against the original
head(ewR$eyewidth.ss)
head(ew$sizeStandardize) # these are very close but slightly different. I think that's because Chelsea said she rounded the commonWithinGroupSlope parameter to 0.62, whereas I'm using the whole thing. Indeed, I tried it out with the parameter rounded to 0.62, and then they're the same out to like 3 decimals.

# At this point, Chelsea also created a model for lake-specific size-standardized values. But she says that those values were not used in future calculations, so I'm going to skip doing that for now.

# I think that this file is now ready to write out and use in other analyses.

# Write out the data ------------------------------------------------------
write.csv(ewR, file = here("data", "outputs", "eyewidthsFINAL.csv"), row.names = F)
