# Script to create Identifiers_Update_2020.txt and Replicates Identifiers.txt, which are used in analysis.R
# Created by Kaija Gahm on 19 February 2021, updated 28 June 2021

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
  table() # all got measured 5 times. The first measurement is used for analysis, and then the four other measurements are used to get measurement error. Because we used the first measurement for our analyses, let's extract replicates 2:5 for the measurement error.

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

# Initialize Identifiers_Update_2020 --------------------------------------
iuR <- data.frame(imageID = filenames)

# Initialize Replicates Identifiers --------------------------------------
riR <- replicates %>%
  rename("imageID" = imageFile,
         "replicateNo" = replic)

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

# Write out Identifiers_Update_2020.txt -----------------------------------
write.table(iuR, here("data", "outputs", "Identifiers_Update_2020.txt"), 
            sep = ",", row.names = F)

# Write out Replicates Identifiers.txt ------------------------------------
write.table(riR, here("data", "outputs", "ReplicatesIdentifiers.txt"), 
            sep = ",", row.names = F)