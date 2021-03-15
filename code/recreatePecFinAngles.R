# Script to recreate pec fin angle data
# Created by Kaija Gahm on 10 March 2021

# Load packages -----------------------------------------------------------
library(dplyr) # for data wrangling
library(tidyr) # for data wrangling
library(stringr) # for text string manipulation
library(here) # for file paths
source(here("code", "dbUtil.R")) # for the dbTable() function
library(RSQLite) # for db connection

# Connect to the database -------------------------------------------------
dbdir <- here("data") # directory where the db is stored
db <- "MFEdb_20210305.db" # name of db

# Load the original pec fin angle data for comparison ---------------------
pfa <- read.csv(here("data", "unclassified", "Pec Fin Angles_ORIGINAL.csv"))

# Grab FISH_MORPHOMETRICS -------------------------------------------------
fm <- dbTable("fish_morphometrics")

# Grab lakeInfo -----------------------------------------------------------
lakeInfo <- read.csv(here("data", "outputs", "lakeInfo_wBins.csv"))

# Initialize recreated df -------------------------------------------------
pfaR <- fm %>%
  select(imageFile, parameter, parameterValue) %>%
  rename("fishID" = imageFile) %>%
  filter(parameter %in% c("X13", "Y13", "X14", "Y14", "pecFinInsertionAngle")) %>%
  group_by(fishID, parameter) %>%
  slice(1) %>% # take the first measurement when there were multiple measurements.
  pivot_wider(id_cols = fishID, names_from = "parameter", values_from = "parameterValue")

# Are all the fish represented?
all(pfa$fishID %in% pfaR$fishID) # good!

# Limit it to the fish contained in the original file
pfaR <- pfaR %>%
  filter(fishID %in% pfa$fishID)

nrow(pfaR) == nrow(pfa)
all(pfaR$fishID == pfa$fishID)
round(pfaR$X13, 2) == round(pfa$X13, 2) # still many not equal.
tail(pfa$X13)
tail(pfaR$X13) # hm... Could this be another case of the multiple-fish problem?
fm %>% filter(imageFile == "TW FN 015.jpg", parameter == "X13") # huh, no, there's only one X13 for e.g. this one.
fm %>% filter(imageFile == "TW FN 014.jpg", parameter == "X13") # only one for this one too. So why are these different?

# XXX START HERE

# lakeID ------------------------------------------------------------------
pfaR <- pfaR %>%
  mutate(lakeID = factor(str_extract(fishID, "^[A-Z]+(?=\\s)"))) %>%
  mutate(lakeID = forcats::fct_recode(lakeID,
                                      "Papoose" = "PP",
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
                                      "Papoose" = "PS",
                                      "Red_Bass" = "RB",
                                      "Red_Bass" = "RS",
                                      "Squaw" = "SQ",
                                      "Towanda" = "TO",
                                      "Towanda" = "TW"))

table(pfaR$lakeID, exclude = NULL) # looks good.
table(pfa$lakeID, exclude = NULL) # good, the counts line up and there are no NA's. I'm using lake names to avoid incorrect abbreviations.

# Add DOC and basin ------------------------------------------------------
pfaR <- pfaR %>%
  left_join(lakeInfo %>%
              select(lakeName, DOC, basin),
            by = c("lakeID" = "lakeName"))

## Check that the DOC values went through
all(pfaR$DOC == pfa$DOC) # yay!

# Calculate fish standard length ------------------------------------------
# Because the photo names are different for the pec fins vs. for the full fish bodies, I need to join through FISH_MORPHOMETRICS. We need a three-way lookup table...

## Get pec fin image file names
pecFinPhotos <- pfaR$fishID

## Get fishID's for those image files
fishIDs <- fm %>%
  filter(parameter == "pecFinInsertionAngle",
         imageFile %in% pecFinPhotos) %>%
  select(fishID, imageFile) %>%
  rename("pfaImageFile" = imageFile) %>%
  distinct()

## Now use those same fishID's to grab the bodyImageFiles and the parameters we need
params <- fm %>%
  filter(parameter %in% c("X7", "Y7", "X1", "Y1"),
         fishID %in% fishIDs$fishID) %>%
  select(fishID, imageFile, parameter, parameterValue) %>%
  group_by(imageFile, parameter) %>%
  slice(1) %>%
  ungroup() %>%
  pivot_wider(names_from = parameter, values_from = parameterValue) %>%
  rename("bodyImageFile" = imageFile)

nrow(params) == nrow(fishIDs)
nrow(params) # 417

## Join the two together
fsl <- left_join(fishIDs, params, by = "fishID")
nrow(fsl) # still 417, good!

## Compute lengths
fsl <- fsl %>%
  mutate(fishStdLength = sqrt((X7-X1)^2+(Y7-Y1)^2)) %>%
  select(pfaImageFile, fishStdLength) %>%
  rename("fishID" = pfaImageFile)

all(fsl$fishID %in% pfaR$fishID) # now we have body lengths for all of the fish! Yay!

## join
pfaR <- pfaR %>%
  left_join(fsl, by = "fishID")

## Check
pfaR$fishStdLength %>% head(10)
pfa$fishStdLength %>% head(10)
# Ok good, these are really really similar.

# Maybe we don't need the X and Y parameters?
head(pfa$ang_deg)
head(pfaR$pecFinInsertionAngle) 
all(pfa$ang_deg == pfaR$pecFinInsertionAngle) # ok yeah these are identical. So presumably the angle measure was taken from this document and put into notoriousBLG. But if that's the case, it was presumably calculated from the X and Y landmarks, so I would like to be able to recreate that. XXX ASK CHELSEA!

# Check to see if relationship with size to see if the angle needs to be size-standardized
ggplot(pfaR, aes(x = fishStdLength, 
                 y = pecFinInsertionAngle)) + 
  geom_point(aes(colour = lakeID)) 

# No size relationship; can leave it as is. No size-standardization needed.

# Write out the data ------------------------------------------------------
write.csv(pfaR, file = here("data", "outputs", "PecFinAnglesFINAL.csv"), row.names = F)
