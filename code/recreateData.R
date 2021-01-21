# Script to create intermediate data files
# Created by Kaija Gahm in January 2021

# Load packages -----------------------------------------------------------
library(RSQLite) # for database connections
library(dplyr) # for data wrangling
library(tidyr) # for data wrangling
library(shadow) # degree/radian conversions
source(here("code", "dbUtil.R")) # for the dbTable() function

# Connect to the database -------------------------------------------------
dbdir <- here("data") # directory where the db is stored
db <- "MFEdb_20201125.db" # name of db

# Pull in FISH_MORPHOMETRICS ----------------------------------------------
fm <- dbTable("fish_morphometrics")

# Pull in lake DOC data ---------------------------------------------------
l <- read.csv(here("data", "outputs", "lakeInfo_wBins.csv"))

# Pec Fin Angles ----------------------------------------------------------
# Based on looking at the raw Pec Fin Angles.csv, it looks like Chelsea took the first measurement for any fish that were measured more than once (i.e. the "replicates" fish)
pfa <- fm %>%
  filter(parameter %in% c("X13", "Y13", "X14", "Y14", "X7", "Y7", "X1", "Y1")) %>%
  select(imageFile, parameter, parameterValue) %>%
  rename("fishID" = imageFile) %>%
  group_by(fishID, parameter) %>%
  slice(1) %>%
  ungroup() %>%
  tidyr::pivot_wider(., names_from = "parameter", 
              values_from = "parameterValue") %>%
  mutate(fishStdLength = sqrt((X7-X1)^2 + (Y7-Y1)^2), # This isn't coming out exactly the same. Check number of digits when Chelsea sends me the excel file.
         lakeID = substr(fishID, 1, 2),
         captureMethod = stringr::str_extract(fishID, "FN|AN|Electro") %>%
           forcats::fct_recode(., 
                      "Angling" = "AN",
                      "Fyke_Net" = "FN",
                      "Electrofishing" = "Electro")) %>% 
  left_join(l %>% select(lakeID, basin, DOC, DOClevel),
            by = "lakeID") %>% 
  select(-c("X1", "Y1", "X7", "Y7")) %>%
  mutate(X = X14-X13,
         Y = Y14-Y13,
         hyp_dist = sqrt(X^2 + Y^2),
         opp_dist = abs(X),
         sinAngle = opp_dist/hyp_dist,
         ang_deg = shadow::rad2deg(asin(sinAngle))) %>%
  as.data.frame() 

