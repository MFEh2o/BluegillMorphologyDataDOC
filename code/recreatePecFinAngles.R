# Script to recreate pectoral fin insertion angle data
# Created by Kaija Gahm on 10 March 2021

# Load packages -----------------------------------------------------------
library(dplyr) # for data wrangling
library(tidyr) # for data wrangling
library(stringr) # for text string manipulation
library(here) # for file paths
source(here("code", "dbUtil.R")) # for the dbTable() function
source(here("code", "defs.R"))
library(RSQLite) # for db connection

# Connect to the database -------------------------------------------------
dbdir <- here("data") # directory where the db is stored
db <- "MFEdb_20210423.db" # name of db

# Grab FISH_MORPHOMETRICS -------------------------------------------------
fm <- dbTable("fish_morphometrics")

# Grab lakeInfo -----------------------------------------------------------
lakeInfo <- read.csv(here("data", "outputs", "Lake_Info_2020wBasins.csv"))

# Initialize recreated df -------------------------------------------------
pfaR <- fm %>%
  select(lakeID, imageFile, parameter, parameterValue) %>%
  rename("fishID" = imageFile) %>%
  filter(parameter %in% c("pecFinInsertionAngle", "stdLength")) %>%
  group_by(lakeID, fishID, parameter) %>%
  slice(1) %>% # take the first measurement when the fish was measured more than once.
  pivot_wider(id_cols = c("lakeID", "fishID"), names_from = "parameter", values_from = "parameterValue")

# Add DOC and basin ------------------------------------------------------
pfaR <- pfaR %>%
  left_join(lakeInfo %>%
              select(lakeID, DOC, basin),
            by = "lakeID")

# Check for size relationship between angle and fish length ---------------
# Does the angle need to be size-standardized?
ggplot(pfaR, aes(x = stdLength, 
                 y = pecFinInsertionAngle)) + 
  geom_point(aes(col = lakeID)) 

# No size relationship; can leave it as is. No size-standardization needed.

# Write out the data ------------------------------------------------------
write.csv(pfaR, file = here("data", "outputs", "PecFinAnglesFINAL.csv"), row.names = F)
