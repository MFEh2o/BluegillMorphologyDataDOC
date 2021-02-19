# Script to create Identifiers_Update_2020.txt
# Created by Kaija Gahm on 19 February 2021

# Load packages -----------------------------------------------------------
library(RSQLite) # for database connections
library(dplyr) # for data wrangling
library(tidyr) # for data wrangling
library(stringr) # for text string manipulation
library(here) # for file paths
source(here("code", "dbUtil.R")) # for the dbTable() function
library(openxlsx)

# Load the original Identifiers_Update_2020.txt for comparison ------------
iu <- read.table(here("data", "dontNeed", "Identifiers_Update_2020_ORIGINAL.txt"), sep = "\t", header = T)

## Remove the underscores in the file names
iu <- iu %>%
  mutate(imageID = str_replace_all(imageID, "_", " "))

# Load other data ---------------------------------------------------------
filenames <- read.table(here("data", "inputs", "fishBodyPhotos_fileNames.txt"), sep = "\t", header = F) %>%
  filter(V1 != "fishBodyPhotos_fileNames.txt") %>%
  pull(V1)

lakeInfo <- read.csv(here("data", "outputs", "lakeInfo_wBins.csv"), header = T)

# Check that these match up with the ones in iu
all(filenames %in% iu$imageID) & all(iu$imageID %in% filenames) # awesome!
all(filenames == iu$imageID) # they're in the same order, too!

# Initialize recreated version of iu --------------------------------------
iuR <- data.frame(imageID = filenames)

# Create the lakeID column ------------------------------------------------
iuR <- iuR %>%
  mutate(lakeID = factor(str_extract(imageID, "^[A-Z]+(?=\\s)"))) %>%
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

## Check that that went through
all(iuR$lakeID == iu$lakeID) # good!

# Add DOC and basin ------------------------------------------------------
iuR <- iuR %>%
  left_join(lakeInfo %>%
              select(lakeName, "lakeDOC" = DOC, 
                     "DOC" = DOClevel, 
                     "DOCrange" = DOCbin,
                     "Basin" = basin),
            by = c("lakeID" = "lakeName"))

## Check that the DOC values went through
all(iuR$lakeDOC == iu$lakeDOC)
iuR %>% filter(lakeDOC != iu$lakeDOC)
data.frame(iuR$lakeDOC, iu$lakeDOC) %>%
  distinct() # ok, this ends up being the same, just with different precision. But I don't think we actually really use the DOC values from the identifiers sheet anyway, so it shouldn't be a problem.

## Check the ranges
iuR %>%
  select(lakeID, lakeDOC, DOC, DOCrange) %>%
  distinct()

iu %>%
  select(lakeID, lakeDOC, DOC, DOCrange) %>%
  distinct() # looks good!

# Add capture method ------------------------------------------------------
iuR <- iuR %>%
  mutate(captureMethod = str_extract(imageID, "AN|FN|Electro")) %>%
  mutate(captureMethod = forcats::fct_recode(captureMethod,
                                    "Angling" = "AN",
                                    "Electrofishing" = "Electro",
                                    "Fyke_Net" = "FN"))

table(iuR$captureMethod, exclude = NULL)
table(iu$captureMethod, exclude = NULL) # good, these match.

# Write out Identifiers_Update_2020.txt -----------------------------------
write.table(iuR, here("data", "outputs", "Identifiers_Update_2020.txt"))