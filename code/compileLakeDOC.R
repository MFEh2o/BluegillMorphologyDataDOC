# Script to assign DOC values for each lake
# Created by Kaija Gahm on 15 April 2021
# Updated 28 April 2021

# Rules for how we chose DOC values:
# 1. If the DOC value from a lake is available in the Craig et al. 2017 paper, use that. Those values were computed from several samples and have methodology available.
# 2. If not #1, but the DOC is available from NTL-LTER Biocomplexity data, then use that value.
# 3. If not #2, then use the value available in the MFE database.

# The lakes we're using are: BA (Bay), BH (Birch), CR (Crampton), FD (Found), HB (Hummingbird), LC (Little Crooked), LT (Lost), MC (McCullough), MS (Muskellunge), OB (Oxbow), PS (Papoose), RS (Red Bass), SQ (Squaw), TO (Towanda).
lakeIDs <- c("BA", "BH", "CR", "FD", "HB", "LC", "LT", "MC", "MS", "OB", "PS", "RS", "SQ", "TO")

# Load packages -----------------------------------------------------------
library(tidyverse) # for data wrangling
library(here) # for file paths
library(RSQLite) # for database connection

# Connect to the database -------------------------------------------------
con <- dbConnect(SQLite(), here("data", "MFEdb_20210423.db"))

# Get lake data from the MFE database -------------------------------------
l <- dbReadTable(con, "lakes")

# Grab the lakes we want --------------------------------------------------
lakeInfo <- l %>%
  filter(lakeID %in% lakeIDs) %>%
  select(lakeID, lakeName, lat, long, WBIC)

# Add DOC for the ones that are in the Craig paper ------------------------
lakeInfo <- lakeInfo %>%
  mutate(DOC = case_when(lakeID == "CR" ~ 5.0,
                         lakeID == "BA" ~ 7.4,
                         lakeID == "BH" ~ 12.5,
                         lakeID == "MC" ~ 14.3,
                         lakeID == "RS" ~ 18.9,
                         lakeID == "HB" ~ 24.5),
         DOCSource = case_when(!is.na(DOC) ~ "Craig et al"))

# Load NTL Biocomplexity data ---------------------------------------------
biocomp <- read.csv(here("data", "inputs", "ntl41_v1_0.csv"))

# The lakeID's here are not necessarily the same as ours.
unique(biocomp$lakename) %>% sort()

lakeNames <- c("Big Muskellunge Lake", "Birch Lake", "Found Lake", "Little Crooked Lake", "Lost Lake", "Muskellunge Lake", "Oxbow Lake", "Papoose Lake", "Towanda Lake")

# Pull just the data we need, and change the lakeID's to match our database
biocomp <- biocomp %>%
  filter(lakename %in% lakeNames) %>%
  select("lakeID" = lakeid, "lakeName" = lakename, 
         "DOC" = doc, flag_doc, sampledate) %>%
  filter(DOC != "NA") %>%
  mutate(lakeID = case_when(lakeName == "Found Lake" ~ "FD",
                            lakeName == "Lost Lake" ~ "LT",
                            lakeName == "Towanda Lake" ~ "TO",
                            lakeName == "Birch Lake" ~ "BH",
                            lakeName == "Papoose Lake" ~ "PS",
                            lakeName == "Muskellunge Lake" ~ "MS",
                            lakeName == "Oxbow Lake" ~ "OB",
                            TRUE ~ lakeID))

all(biocomp$lakeID %in% lakeIDs) # TRUE

# Get the average DOC value for each lake
toJoin <- biocomp %>%
  group_by(lakeID) %>%
  summarize(DOC = mean(DOC)) %>%
  mutate(DOCSource = "NTL biocomplexity data")

# Join those values to the other lake info
lakeInfo <- lakeInfo %>%
  rows_update(x = .,
              y = toJoin,
              by = "lakeID")

# Now only Squaw lake is left. For that one, let's take the value from the MFE database.
wc <- tbl(con, "water_chem") %>%
  filter(parameter == "DOC",
         lakeID == "SQ") %>%
  select(sampleID, dateSample, parameterValue, flag) %>%
  collect()
wc # There are only two values for SQ. Let's use the average of the two. 

wc <- wc %>%
  summarize(DOC = mean(as.numeric(parameterValue))/1000) %>%
  pull(DOC)

lakeInfo[lakeInfo$lakeID == "SQ", c("DOC", "DOCSource")] <- c(wc, "MFE database")

# Write out the lake data
write.csv(lakeInfo, here("data", "outputs", "lakeInfo.csv"))

