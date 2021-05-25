# Script to create Identifiers_Update_2020.txt and Replicates Identifiers.txt
# Created by Kaija Gahm on 19 February 2021, updated 20 May 2021

# Load packages -----------------------------------------------------------
library(dplyr) # for data wrangling
library(tidyr) # for data wrangling
library(stringr) # for text string manipulation
library(here) # for file paths
library(openxlsx)
library(RSQLite) # for db connection
source(here("code", "dbUtil.R"))

# Connect to the database -------------------------------------------------
dbdir <- here("data") # directory where the db is stored
db <- "MFEdb_20210423.db" # name of db

# Load the original text files for comparison ------------
iu <- read.table(here("archived", "data", "Identifiers_Update_2020_ORIGINAL.txt"), sep = "\t", header = T)
ri <- read.table(here("archived", "data", "Replicates Identifiers_ORIGINAL.txt"), sep = "\t", header = T)

## Remove the underscores in the file names for iu
iu <- iu %>%
  mutate(imageID = str_replace_all(imageID, "_", " "))

# Load other data ---------------------------------------------------------
filenames <- read.table(here("data", "inputs", "fishBodyPhotos_fileNames.txt"), sep = "\t", header = F) %>%
  filter(V1 != "fishBodyPhotos_fileNames.txt") %>%
  pull(V1)

lakeInfo <- read.csv(here("data", "outputs", "Lake_Info_2020wBasins.csv"), header = T)

# get the landmark data from the database
fm <- dbTable("fish_morphometrics")

# grab only the fish where replicates were taken
replicates <- fm %>% 
  filter(parameter == "X1") %>%
  group_by(imageFile) %>%
  filter(n() > 1) %>% # only the ones that had replicates taken
  arrange(imageFile) %>% # sort in order
  ungroup()

# How many times did each fish get measured?
replicates %>%
  group_by(imageFile) %>%
  summarize(n = n()) %>%
  pull(n) %>%
  table() # all got measured 5 times. That's weird, I thought it was 4. 

# I guess maybe they were measured one time (the one that got used for future analyses) and then four *more* times to get measurement error. Because we used the first measurement for our analyses, let's extract replicates 2:5 for the measurement error.

replicates <- replicates %>%
  group_by(imageFile) %>%
  slice(2:5)

nrow(replicates) # 240
replicates$imageFile # these are ordered with all the replicates together. Let's order them to match ri, so that we go through all image files for one replicate before going through the next one.

replicates <- replicates %>%
  group_by(imageFile) %>%
  mutate(replic = 1:n()) %>%
  ungroup() %>%
  select(imageFile, replic) %>%
  arrange(replic, imageFile)

all(replicates$imageFile == ri$imageID) #okay, great, these are in the right order for recreating ri.

# Check that these match up with the ones in iu
all(filenames %in% iu$imageID) & all(iu$imageID %in% filenames) # awesome!
all(ri$imageID %in% filenames) # not the other way around, but that's okay
all(filenames == iu$imageID) # they're in the same order, too, for iu!

# Initialize recreated version of iu --------------------------------------
iuR <- data.frame(imageID = filenames)

# Initialize recreated version of ri --------------------------------------
riR <- replicates %>%
  rename("imageID" = imageFile,
         "replicateNo" = replic)

all(riR$imageID == ri$imageID) # all match, yay!
all(riR$replicateNo == ri$replicateNo) # all match, yay!

# Create the lakeID column ------------------------------------------------
iuR <- iuR %>%
  mutate(lakeID = word(imageID, 1, 1, sep = "\\s")) %>%
  mutate(lakeID = forcats::fct_recode(lakeID,
                                      "MS" = "MK",
                                      "PS" = "PP",
                                      "RS" = "RB",
                                      "TO" = "TW"))
table(iuR$lakeID) # looks good

# Add DOC and basin ------------------------------------------------------
iuR <- iuR %>%
  left_join(lakeInfo %>%
              select(lakeID, "lakeDOC" = DOC, 
                     "Basin" = basin),
            by = "lakeID")

table(iuR$lakeDOC, exclude = NULL)

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

# Write out Replicates Identifiers.txt ------------------------------------
write.table(riR, here("data", "outputs", "Replicates Identifiers.txt"), sep = ",", row.names = F)

