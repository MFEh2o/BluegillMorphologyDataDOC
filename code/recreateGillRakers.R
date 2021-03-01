# Script to recreate Gill_Rakers_2018_Final.csv
# Created by Kaija Gahm on 25 February 2021

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

# Load the original Gill_Rakers_2018_Final.csv for comparison -------------
gr <- read.csv(here("data", "unclassified", "Gill_Rakers_2018_Final_ORIGINAL.csv")) 

# Grab FISH_MORPHOMETRICS -------------------------------------------------
fm <- dbTable("fish_morphometrics")

# Grab lakeInfo -----------------------------------------------------------
lakeInfo <- read.csv(here("data", "outputs", "lakeInfo_wBins.csv"))

# Initialize recreated df -------------------------------------------------
grR <- fm %>%
  select(imageFile, parameter, parameterValue) %>%
  rename("fishID" = imageFile)

## filter to include only the gill raker parameters, and recode them to match the column names in the original Gill_Rakers_2018 sheet
grR <- grR %>%
  filter(grepl("Raker", parameter)) %>%
  mutate(parameter = str_replace(parameter, "length", "length_"),
         parameter = str_replace(parameter, "spaceRaker", "space")) %>%
  pivot_wider(id_cols = "fishID",
              names_from = "parameter",
              values_from = "parameterValue") %>%
  rename("total_RakerNum" = "totalRakerCount")

all(names(grR) %in% names(gr)) # yes, all names show up.

# Are all the fish represented?
all(gr$fishID %in% grR$fishID) # good!

# Limit it to the fish contained in the original file
grR <- grR %>%
  filter(fishID %in% gr$fishID)

# check a few columns and dimensions to make sure they match
nrow(grR) == nrow(gr)
all(grR$fishID == gr$fishID)
all(grR$length_Raker1 == gr$length_Raker1)
all(grR$space1 == gr$space1) # all good!

# lakeID ------------------------------------------------------------------
grR <- grR %>%
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

table(grR$lakeID, exclude = NULL)
table(gr$lakeID, exclude = NULL) # good, the counts line up and there are no NA's. I'm using lake names to avoid incorrect abbreviations.

# Add DOC and basin ------------------------------------------------------
grR <- grR %>%
  left_join(lakeInfo %>%
              select(lakeName, "lakeDOC" = DOC, DOClevel, basin, DOCbin,
                     "area" = surface_Area.ha., 
                     "max_Depth" = max_Depth.m.),
            by = c("lakeID" = "lakeName"))

## Check that the DOC values went through
all(grR$lakeDOC == gr$lakeDOC) # yay!

# Not adding capture method because it doesn't get used in analysis.

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
  filter(fishID %in% grR$fishID)

## join
grR <- grR %>%
  left_join(fsl %>%
              rename("fishSL" = "fishStdLength"), 
            by = "fishID")

## Check
grR$fishSL %>% head(10)
gr$fishSL %>% head(10)# very close! Good enough.

# Compute average raker lengths and spaces --------------------------------
grR <- grR %>%
  rowwise() %>%
  mutate(avgRakerLength = mean(c_across(contains("length_Raker"))),
         avgRakerSpace = mean(c_across(contains("space"))),
         avgL_4.7 = mean(c_across(c(length_Raker4, length_Raker5, length_Raker6, length_Raker7))),
         avgS_4.6 = mean(c_across(c(space4, space5, space6))))

## check that these averages look good
head(grR$avgL_4.7)
head(gr$avgL_4.7) # okay, great! These look good. I also checked all of them, and some are off by slight amounts, but it's definitely just a rounding thing.
# KG: double check that this is actually just a rounding error. 

# Size-standardize the averages -------------------------------------------
# Looking at the Review April 2020 script, it seems that the only size-standardized measurements that are actually used in the model are the avgL_4.7 measurement and the avgS_4.6 measurement, both size-standardized independent of lakeID (so, size-standardized using the full data set, as seen in recreateEyewidths.R).
# So, I'm only going to perform size-standardizations for those columns.

# Average Raker Length (Rakers 4-7)
## To get the common within-group slope to be used in size-standardization, we need to exclude an interaction effect with lakeID. So we'll fit this model:
size.rem.avgL <- lm(log10(avgL_4.7) ~ log10(fishSL)+lakeID, data = grR)
summary(size.rem.avgL) # now we can extract the common within-group slope, which is the `Estimate` parameter for `log10(fishSL)`: 0.398793

commonWithinGroupSlope <- coef(size.rem.avgL)[2] %>% unname() # programmatically grab that coefficient
meanFishSize <- mean(grR$fishSL) # compute the mean `fishSL` across all fish in the data set

## Now add a column to grR for size-standardized avgL_4.7
grR <- grR %>%
  mutate(avgL2_ss = avgL_4.7*(meanFishSize/fishSL)^commonWithinGroupSlope)

## check it against the original
head(grR$avgL2_ss)
head(gr$avgL2_ss) # These values are slightly different than the values in Chelsea's csv. She looked into it, and explains that she accidentally used coef(size.rem.avgL)[1] (the intercept value) instead of coef(size.rem.avgL)[2] (the log10(fishSL) value) as the commonWithinGroupSlope value Since the former is 0.424219 and the latter is 0.398793, the resulting values are only slightly different. The same mistake was made for the raker spaces (see below), resulting in a much greater difference in values since the intercept and slope values were much more different in that model.

# Average Raker Space (Spaces 4-6)
## To get the common within-group slope to be used in size-standardization, we need to exclude an interaction effect with lakeID. So we'll fit this model:
size.rem.avgS <- lm(log10(avgS_4.6) ~ log10(fishSL)+lakeID, data = grR)
summary(size.rem.avgS) # now we can extract the common within-group slope, which is the `Estimate` parameter for `log10(fishSL)`: 0.72846

commonWithinGroupSlope <- coef(size.rem.avgS)[2] %>% unname() # programmatically grab that coefficient
# We don't need to re-compute meanFishSize, since it'll just be the same as above.

## Now add a column to grR for size-standardized avgS_4.6
grR <- grR %>%
  mutate(avgS2_ss = avgS_4.6*(meanFishSize/fishSL)^commonWithinGroupSlope)

## check it against the original
head(grR$avgS2_ss)
head(gr$avgS2_ss) # These values are extremely different. See the explanation above for raker lengths: in her original work, Chelsea took the wrong coefficient, the intercept instead of the slope. Since the former, for this model, is -0.59703633 and the latter is 0.72846002, the values are understandably hugely different.

# Check names between grR and gr.
names(gr)[!(names(gr) %in% names(grR))]
# Okay, we still have not recreated all the column names, but that's okay. capture_Method isn't necessary for analyses. SS.Length, SS.Space, and SS.Count were also not used in subsequent analyses, so I haven't performed size-standardizations on those. lakeSlope, avgrakerlengthSS_bylake, avgrakerspaceSS_bylake, and rakercountSS_bylake are all by-lake size-standardized values, but Chelsea says that these were part of some exploratory analyses and weren't used in subsequent models. fitted_L, fitted_S, and fitted_C will be added in the 'Review April 2020_KGedit.R' script. So we're ready to write out the data.

# Write out the data ------------------------------------------------------
write.csv(grR, file = here("data", "outputs", "Gill_Rakers_2018_Final.csv"), row.names = F)