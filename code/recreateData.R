# Script to create intermediate data files
# Created by Kaija Gahm in January 2021

# Load packages -----------------------------------------------------------
library(RSQLite) # for database connections
library(dplyr) # for data wrangling
library(tidyr) # for data wrangling
library(stringr) # for text string manipulation
library(shadow) # degree/radian conversions
library(here) # for file paths
source(here("code", "dbUtil.R")) # for the dbTable() function
library(kgUtils)

# Connect to the database -------------------------------------------------
dbdir <- here("data") # directory where the db is stored
db <- "MFEdb_20201125.db" # name of db

# Pull in FISH_MORPHOMETRICS ----------------------------------------------
fm <- dbTable("fish_morphometrics")
head(fm)

# Pull in identifiers -----------------------------------------------------
## Want to compare this to the data in FISH_MORPHOMETRICS
iu <- read.table(here("data", "unclassified", "Identifiers_Update_2020.txt"), 
                 sep = "\t", header = T)
nrow(iu) #417 rows
length(unique(iu$imageID)) # 417 unique images

# Going to remove the underscores from iu image identifiers so we can try to reconcile them
iu <- iu %>%
  mutate(imageID = str_replace_all(imageID, "_", " "))

# Compare the imageID's in iu to the imageID's in FISH_MORPHOMETRICS (fm)
length(unique(iu$imageID)) # still 417, good
length(unique(fm$imageFile)) # 635. What's going on?

fm %>%
  filter(!(imageFile %in% iu$imageID)) %>%
  pull(imageFile) %>%
  unique() # okay, so it looks like the ones that aren't in IU are the ones where the order of the date and the photo number is reversed. Can we un-reverse them? Are those the same fish?

# Let's just see what happens if we do
lookup <- fm %>%
  select(imageFile) %>%
  distinct()

## Separate lookup into problems and not problems
problems <- lookup %>%
  filter(!(imageFile %in% iu$imageID))

notProblems <- lookup %>%
  filter(imageFile %in% iu$imageID)

## For problems, switch the order of the parts
problems <- problems %>%
  mutate(date = str_extract(imageFile, "\\d{1,2}\\-\\d{1,2}"),
         number = str_extract(imageFile, "\\d{3}(?=\\.jpg)"),
         prefix = str_replace(imageFile, pattern = date, replacement = ""),
         prefix = str_replace(prefix, pattern = number, replacement = ""),
         prefix = str_replace(prefix, pattern = "\\s\\s\\.jpg", replacement = ""),
         fixed = paste0(prefix, " ", number, " ", date, ".jpg") # have to fix the NA that gets introduced here.
  )

## Fix any that include NA's
problems %>%
  filter(grepl("NA", fixed))

problems$fixed[problems$imageFile == "RS FN4 002.jpg"] <- "RB FN4 002.jpg" # Red Bass is supposed to be "RS", but it seems like it was encoded in IU as "RB". RB is actually Raspberry, but in the case of this project, it's supposed to refer to "Red Bass". I will have to make this change across the board in FM and in all the other sheets. *Definitely* in FM, since that's in the database.

fm <- fm %>%
  left_join(problems %>% select(imageFile, fixed),
            by = "imageFile")

fm <- fm %>%
  mutate(imageFile = case_when(!is.na(fixed) ~ fixed,
                               TRUE ~ imageFile))

## Now, how many of these match up with IU?
unique(fm$imageFile) %>% length()  # now we're down to 581 instead of 635. Still concerning.

fm %>%
  select(imageFile) %>%
  distinct() %>%
  filter(!(imageFile %in% iu$imageFile)) # a whole bunch of these are still not represented in iu.



# Pull in lake DOC data ---------------------------------------------------
l <- read.csv(here("data", "outputs", "lakeInfo_wBins.csv"))


# Fish standard lengths and basic info ------------------------------------
fsl <- fm %>%
  filter(parameter %in% c("X7", "Y7", "X1", "Y1")) %>%
  select(imageFile, parameter, parameterValue) %>%
  group_by(imageFile, parameter) %>% # because standard length was calculated based on just the first measurement, we group by imageFile and parameter and select the first row.
  slice(1) %>%
  ungroup() %>%
  pivot_wider(names_from = parameter, values_from = parameterValue) %>%
  rename("fishID" = imageFile) %>%
  mutate(fishStdLength = sqrt((X7-X1)^2+(Y7-Y1)^2),
         lakeID = substr(fishID, 1, 2),
         captureMethod = stringr::str_extract(fishID, "AN|FN|Electro") %>%
           forcats::fct_recode(., 
                               "Angling" = "AN",
                               "Fyke_Net" = "FN",
                               "Electrofishing" = "Electro")) %>%
  select(-c("X1", "X7", "Y1", "Y7")) %>%
  mutate(lakeID = forcats::fct_recode(lakeID, # some of the lakeID's were wrong, so we will fix them.
                                      "PS" = "PP",
                                      "MS" = "MK"))
  
# Check standard lengths against the ones from the gill rakers spreadsheet
gro <- read.csv(here("data", "unclassified", 
                     "Gill_Rakers_2018_Final_ORIGINAL.csv")) %>%
  select(fishID, fishSL)
head(gro)
head(fsl) # looks like we're accurate out to the second or third decimal place. 

all(gro$fishID %in% fsl$fishID) # good.
all(fsl$fishID %in% gro$fishID) # uh oh.
fsl %>% # sure enough, it looks like there are a few fish in fsl that don't show up in gro. Why?
  filter(!(fishID %in% gro$fishID)) %>% 
  pull(fishID)

# - check digits and rounding for the fishStdLength measurements, since they came out slightly different than Chelsea's (Chris thinks this is not a concern, since it's probably just rounding error, and they are accurate to a couple decimal places)

# Pec Fin Angles ----------------------------------------------------------
# Based on looking at the raw Pec Fin Angles.csv, it looks like Chelsea took only the first measurement for fish that were measured more than once (i.e. the "replicates" fish)
pfa <- fm %>%
  filter(parameter %in% c("X13", "Y13", "X14", "Y14")) %>% # landmarks 13 and 14 are on either side of the base of the pectoral fin
  select(imageFile, parameter, parameterValue, replicate) %>% # including replicate allows us to ID the duplicate fish
  rename("fishID" = imageFile) %>%
  tidyr::pivot_wider(id_cols = c("fishID", "replicate"), names_from = "parameter", 
              values_from = "parameterValue") %>%
  left_join(fsl, by = "fishID") %>% # join standard lengths
  left_join(l %>% select(lakeID, basin, DOC, DOClevel), # join lake information
            by = "lakeID") %>% 
  mutate(X = X14-X13, # compute pec fin angle
         Y = Y14-Y13,
         hyp_dist = sqrt(X^2 + Y^2),
         opp_dist = abs(X),
         sinAngle = opp_dist/hyp_dist,
         ang_deg = shadow::rad2deg(asin(sinAngle))) %>%
  as.data.frame() 

pfao <- read.csv(here("data", "unclassified", "Pec Fin Angles_ORIGINAL.csv"))
# ^ use the original to check the new data. Sure enough, all fish are included. 

write.csv(pfa, here("data", "outputs", "Pec Fin Angles.csv"))
## Still to recreate:
# - in order to recreate the `fitted` column, have to figure out what's going on with the various angle columns that are called for in the univariate data script but aren't present in the data. 

# PecFinDataNovember.csv --------------------------------------------------
# According to thesis, pec fin length is the distance between landmarks 14 and 21. (but don't need to calculate this from landmarks; it's already in the data as "pecFinLength")
# pec fin base width is the distance between landmarks 13 and 14.
pfdn <- fm %>%
  filter(parameter %in% c("pecFinLength", "pecFinBaseWidth")) %>%
  select("fishID" = imageFile,
         parameter,
         parameterValue,
         replicate) %>% 
  pivot_wider(id_cols = c("fishID", "replicate"), names_from = "parameter", values_from = "parameterValue") 

# For some reason, the image file/fishID's in PecFinDataNovember are different from the ones in the other sheets: pfdn has date/number, not number/date.

pfdno <- read.csv(here("data", "unclassified", "PecFinDataNovember.csv"))
test <- pfdno %>%
  select(fishID, lakeID, baseWidth, PecFinLengths, fishStdLength) %>% 
  distinct()

pfdn$fishID %in% pfdno$fishID


# SOME INVESTIGATIONS WITH CHELSEA ----------------------------------------
# Read in the photo names so we can ground-truth it against FISH_MORPHOMETRICS. This list of photos is the official list of the final photos that Chelsea used in her analysis.
f <- read.table(here("data", "unclassified", "fishBodyPhotos_fileNames.txt"), sep = ",")
f <- f %>%
  filter(V1 != "fishBodyPhotos_fileNames.txt")
f$V1 %in% iu$imageID

length(unique(f$V1))
f %>%
  group_by(V1) %>%
  filter(n() > 1)
  
# ok, all of the final names are in Identifiers_Update, which is good. And there are no duplicate photo names--we have 417 photos, which is what we should have.

# Let's compare these photos to FISH_MORPHOMETRICS (note that this is after I did the switcharoos of the date vs. number order, above.)
extras <- fm %>%
  filter(!(imageFile %in% f$V1)) %>%
  pull(imageFile) %>%
  unique() %>%
  sort() # sort them alphabetically

# Write these out for Chelsea
write.csv(extras, here("data", "outputs", "extraFishPhotoNames.20210215.csv"))

# Remaining problems:
# The fish for which we have pecFinLength and pecFinBaseWidth measurements are different fish than the ones for which we have landmarks 1 and 7 (for calculating standard length).

# Then have to go and dig around in the univariate analysis script to figure out which model corresponds to which value.

# Gill_Rakers_2018_Final.csv ----------------------------------------------
gr <- fm %>%
  filter(grepl("Raker", parameter)) %>%
  select("fishID" = imageFile,
         parameter, parameterValue, replicate) %>%
  pivot_wider(id_cols = c("fishID", "replicate"),
              names_from = "parameter",
              values_from = "parameterValue") %>%
  rename("total_RakerNum" = totalRakerCount) %>%
  rowwise() %>% 
  mutate(avgRakerLength = mean(c(lengthRaker1, lengthRaker2, lengthRaker3, 
                                 lengthRaker4, lengthRaker5, lengthRaker6, 
                                 lengthRaker7, lengthRaker8, lengthRaker9)),
         avgRakerSpace = mean(c(spaceRaker1, spaceRaker2, spaceRaker3, 
                                spaceRaker4, spaceRaker5, spaceRaker6,
                                spaceRaker7, spaceRaker8))) %>%
  mutate(SS.Length = NA,
         SS.Space = NA,
         SS.Count = NA,
         lakeSlope = NA,
         avgrakerlengthSS_bylake = NA,
         avgrakerspaceSS_bylake = NA,
         rakercountSS_bylake = NA,
         `avgL_4-7` = mean(lengthRaker4, lengthRaker5, lengthRaker6, lengthRaker7),
         `avgS_4-6` = mean(spaceRaker4, spaceRaker5, spaceRaker6),
         avgL2_ss = NA,	
         avgS2_ss = NA)


  




