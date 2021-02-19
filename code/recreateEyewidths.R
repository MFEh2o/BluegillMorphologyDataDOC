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

# lakeID and DOC ----------------------------------------------------------
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
all(ewR$lakeDOC == ew$lakeDOC)
data.frame(iuR$lakeDOC, iu$lakeDOC) %>%
  distinct() # ok, this ends up being the same, just with different precision. But I don't think we actually really use the DOC values from the identifiers sheet anyway, so it shouldn't be a problem.

# add capture method ------------------------------------------------------
ewR <- ewR %>%
  mutate(captureMethod = str_extract(fishID, "AN|FN|Electro")) %>%
  mutate(captureMethod = forcats::fct_recode(captureMethod,
                                             "Angling" = "AN",
                                             "Electrofishing" = "Electro",
                                             "Fyke_Net" = "FN"))

table(ewR$captureMethod, exclude = NULL)
table(ew$captureMethod, exclude = NULL) # good, these match.
all(ewR$captureMethod == ew$captureMethod) # cool.

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

# Now we just need lakeSlope and lakeSS -----------------------------------
# Or, do we actually need these?
# I don't see these at all in the Review April 2020 script, so I am going to leave them out for now.

# Write out the data ------------------------------------------------------
write.csv(ewR, file = here("data", "outputs", "eyewidthsFINAL.csv"), row.names = F)