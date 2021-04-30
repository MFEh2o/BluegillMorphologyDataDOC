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
db <- "MFEdb_20210423.db" # name of db

# Load the original pec fin angle data for comparison ---------------------
pfa <- read.csv(here("data", "unclassified", "Pec Fin Angles_ORIGINAL.csv"))

# Grab FISH_MORPHOMETRICS -------------------------------------------------
fm <- dbTable("fish_morphometrics")

# Grab lakeInfo -----------------------------------------------------------
lakeInfo <- read.csv(here("data", "outputs", "Lake_Info_2020wBasins.csv"))

# Initialize recreated df -------------------------------------------------
pfaR <- fm %>%
  select(imageFile, parameter, parameterValue) %>%
  rename("fishID" = imageFile) %>%
  filter(parameter %in% c("pecFinInsertionAngle", "stdLength")) %>%
  group_by(fishID, parameter) %>%
  slice(1) %>% # take the first measurement when the fish was measured more than once.
  pivot_wider(id_cols = fishID, names_from = "parameter", values_from = "parameterValue")

# Are all the fish represented?
all(pfa$fishID %in% pfaR$fishID) # good!

# Limit it to the fish contained in the original file
pfaR <- pfaR %>%
  filter(fishID %in% pfa$fishID)

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

# Check for size relationship between angle and fish length ---------------
# Does the angle need to be size-standardized?
ggplot(pfaR, aes(x = stdLength, 
                 y = pecFinInsertionAngle)) + 
  geom_point(aes(colour = lakeID)) 

# No size relationship; can leave it as is. No size-standardization needed.

# Write out the data ------------------------------------------------------
write.csv(pfaR, file = here("data", "outputs", "PecFinAnglesFINAL.csv"), row.names = F)
