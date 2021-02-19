# Script to classify lakes into DOC categories
# Written by Kaija Gahm, January 2021

# Load packages -----------------------------------------------------------
library(here)
library(dplyr)

# Load lake info ----------------------------------------------------------
l <- read.csv(here("data", "outputs", "Lake_Info_2020wBasins.csv"))

# Assign DOC categories ----------------------------------------------------
# According to Chelsea, the categories we will use are: 
# 1. For low/high: 10mg/L is the cutoff
# 2. For bins, they're 0-5, 5-10, 10-15, 15-20, 20-25, etc.

l <- l %>%
  mutate(DOClevel = case_when(DOC > 10 ~ "high",
                              TRUE ~ "low"),
         DOCbin = case_when(DOC < 5 ~ "0-5",
                            DOC >= 5 & DOC < 10 ~ "5-10",
                            DOC >= 10 & DOC < 15 ~ "10-15",
                            DOC >= 15 & DOC < 20 ~ "15-20",
                            DOC >= 20 & DOC < 25 ~ "20-25"))

# Write out ---------------------------------------------------------------
write.csv(l, here("data", "outputs", "lakeInfo_wBins.csv"), row.names = F)
