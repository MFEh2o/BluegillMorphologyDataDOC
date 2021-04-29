# Script to create Identifiers_Update_2020.txt
# Created by Kaija Gahm on 19 February 2021

# Load packages -----------------------------------------------------------
library(dplyr) # for data wrangling
library(tidyr) # for data wrangling
library(stringr) # for text string manipulation
library(here) # for file paths
library(openxlsx)

# Load the original Identifiers_Update_2020.txt for comparison ------------
iu <- read.table(here("archived", "data", "Identifiers_Update_2020_ORIGINAL.txt"), sep = "\t", header = T)

## Remove the underscores in the file names
iu <- iu %>%
  mutate(imageID = str_replace_all(imageID, "_", " "))

# Load other data ---------------------------------------------------------
filenames <- read.table(here("data", "inputs", "fishBodyPhotos_fileNames.txt"), sep = "\t", header = F) %>%
  filter(V1 != "fishBodyPhotos_fileNames.txt") %>%
  pull(V1)

lakeInfo <- read.csv(here("data", "outputs", "Lake_Info_2020wBasins.csv"), header = T)

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
                     "Basin" = basin),
            by = c("lakeID" = "lakeName"))

# Add capture method ------------------------------------------------------
iuR <- iuR %>%
  mutate(captureMethod = str_extract(imageID, "AN|FN|Electro")) %>%
  mutate(captureMethod = forcats::fct_recode(captureMethod,
                                    "Angling" = "AN",
                                    "Electrofishing" = "Electro",
                                    "Fyke_Net" = "FN"))

table(iuR$captureMethod, exclude = NULL)
table(iu$captureMethod, exclude = NULL) # good, these match.

# Sex: Note that I talked to Chelsea, and she said that the Sex data is unreliable. There are some discrepancies between the Sex data in Identifiers_Update and the Sex data in notoriousBLG. She says that's because some fish were sexed in the field, while others were sexed in the lab. But there were some mix-ups with death tags that led to that sex information being appended to the wrong fish... Anyway. She says the sex data is not reliable and I should not try to add it to FISH_MORPHOMETRICS. She also didn't end up using it in her analyses. So I have left it off.

# Write out Identifiers_Update_2020.txt -----------------------------------
write.table(iuR, here("data", "outputs", "Identifiers_Update_2020.txt"), sep = ",", row.names = F)
