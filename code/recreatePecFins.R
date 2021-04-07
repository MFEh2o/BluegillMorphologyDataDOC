# Script to recreate pectoral fin data
# Created by Kaija Gahm on 3 March 2021

# Load packages -----------------------------------------------------------
library(dplyr) # for data wrangling
library(tidyr) # for data wrangling
library(stringr) # for text string manipulation
library(here) # for file paths
source(here("code", "dbUtil.R")) # for the dbTable() function
library(RSQLite) # for db connection

# Connect to the database -------------------------------------------------
dbdir <- here("data") # directory where the db is stored
db <- "MFEdb_20210305.db" # name of db

# Load the original pec fin data for comparison ---------------------------
pf <- read.csv(here("data", "unclassified", "PecFinDataNovember_ORIGINAL.csv")) %>%
  select(-imageID)

# Grab FISH_MORPHOMETRICS -------------------------------------------------
fm <- dbTable("fish_morphometrics")

# Grab lakeInfo -----------------------------------------------------------
lakeInfo <- read.csv(here("data", "outputs", "lakeInfo_wBins.csv"))

# Initialize recreated df -------------------------------------------------
pfR <- fm %>%
  select(imageFile, parameter, parameterValue) %>%
  rename("fishID" = imageFile) %>%
  filter(parameter %in% c("pecFinLength", "pecFinBaseWidth")) %>%
  pivot_wider(id_cols = fishID, names_from = "parameter", values_from = "parameterValue")

# Are all the fish represented?
all(pf$fishID %in% pfR$fishID) # good!

# Limit it to the fish contained in the original file
pfR <- pfR %>%
  filter(fishID %in% pf$fishID)

nrow(pfR) == nrow(pf)
all(pfR$fishID == pf$fishID)
all(pfR$pecFinLength == pf$PecFinLengths) # all right, these match!
all(pfR$pecFinBaseWidth == pf$baseWidth) # these also match. Good!

# lakeID ------------------------------------------------------------------
pfR <- pfR %>%
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
                                      "Papoose" = "PS",
                                      "Red_Bass" = "RS",
                                      "Squaw" = "SQ",
                                      "Towanda" = "TO"))

table(pfR$lakeID, exclude = NULL) # looks good.
table(pf$lakeID, exclude = NULL) # good, the counts line up and there are no NA's. I'm using lake names to avoid incorrect abbreviations.

# Add DOC and basin ------------------------------------------------------
pfR <- pfR %>%
  left_join(lakeInfo %>%
              select(lakeName, DOC, 
                     "DOCbin" = DOClevel, basin),
            by = c("lakeID" = "lakeName"))

## Check that the DOC values went through
all(pfR$DOC == pf$DOC) # yay!

# Calculate fish standard length ------------------------------------------
# Because the photo names are different for the pec fins vs. for the full fish bodies, I need to join through FISH_MORPHOMETRICS. We need a three-way lookup table...

## Get pec fin image file names
pecFinPhotos <- pfR$fishID

## Get fishID's for those image files
fishIDs <- fm %>%
  filter(parameter == "pecFinLength",
         imageFile %in% pecFinPhotos) %>%
  select(fishID, imageFile) %>%
  rename("pfImageFile" = imageFile) %>%
  distinct()

## Now use those same fishID's to grab the bodyImageFiles and the parameters we need
params <- fm %>%
  filter(parameter %in% c("X7", "Y7", "X1", "Y1"),
         fishID %in% fishIDs$fishID) %>%
  select(fishID, imageFile, parameter, parameterValue) %>%
  group_by(imageFile, parameter) %>%
  slice(1) %>%
  ungroup() %>%
  pivot_wider(names_from = parameter, values_from = parameterValue) %>%
  rename("bodyImageFile" = imageFile)

nrow(params) == nrow(fishIDs)

## Join the two together
fsl <- left_join(fishIDs, params, by = "fishID")
nrow(fsl) # still 218, good!

## Compute lengths
fsl <- fsl %>%
  mutate(fishStdLength = sqrt((X7-X1)^2+(Y7-Y1)^2)) %>%
  select(-c("fishID")) %>%
  rename("fishID" = pfImageFile)

all(fsl$fishID %in% pfR$fishID) # now we have body lengths for all of the fish! Yay!

## join
pfR <- pfR %>%
  left_join(fsl %>% select(fishID, fishStdLength), by = "fishID")

## Check
pfR$fishStdLength %>% head(10)
pf$fishStdLength %>% head(10)
# These lengths are pretty far off. However, Chelsea and I both scrutinized them in March-April 2021, and Chelsea says she trusts the recomputed values more than the older values (per her 4/7 email). So I'm going to proceed with these recomputed values, even though that will change the model outputs.

# Size-standardize the pec fins -----------------------------------------
#finLengthSS, finBaseSS
## Pec Fin Length
## To get the common within-group slope to be used in size-standardization, we fit a model without an interaction effect with lakeID
length.rem <- lm(log10(pecFinLength) ~ log10(fishStdLength)+lakeID, data = pfR)
summary(length.rem) # now we can extract the common within-group slope, which is the `Estimate` parameter for `log10(fishStdLength)`: 0.771971.

commonWithinGroupSlope <- coef(length.rem)[2] %>% unname() # programmatically grab that coefficient
meanFishSize <- mean(pfR$fishStdLength) # compute the mean `fishStdLength` across all fish in the data set

## Now add a column to pfR for size-standardized eye width
pfR <- pfR %>%
  mutate(finLengthSS = pecFinLength*(meanFishSize/fishStdLength)^commonWithinGroupSlope)

## check it against the original
head(pfR$finLengthSS)
head(pf$finLengthSS) # These are a bit off, but note that the differences come from the discrepancies in fish lengths. As described above, we're proceeding with the re-computed values.

## Pec Fin Base Width
## To get the common within-group slope to be used in size-standardization, we fit a model without an interaction effect with lakeID
width.rem <- lm(log10(pecFinBaseWidth) ~ log10(fishStdLength)+lakeID, data = pfR)
summary(width.rem) # now we can extract the common within-group slope, which is the `Estimate` parameter for `log10(fishStdLength)`: 0.863176.

commonWithinGroupSlope <- coef(width.rem)[2] %>% unname() # programmatically grab that coefficient
# meanFishSize is the same as it was above.

## Now add a column to pfR for size-standardized eye width
pfR <- pfR %>%
  mutate(finBaseSS = pecFinBaseWidth*(meanFishSize/fishStdLength)^commonWithinGroupSlope)

## check it against the original
head(pfR$finBaseSS)
head(pf$finBaseSS) # Once again, there are discrepancies here, with the differences coming from the differences in fish length described above. We're going with the re-computed values.

# At this point in her original scripts, Chelsea also created a model for lake-specific size-standardized values. But she says that those values were not used in future calculations, so I'm going to skip doing that for now.

# Write out the data ------------------------------------------------------
write.csv(pfR, file = here("data", "outputs", "PecFinDataNovemberFINAL.csv"), row.names = F)

# (Below: I wrote out the old and new fish std length values for Chelsea to look at; this helped her determine which values we should trust.)
# # Write out comparison fish sizes for Chelsea to look at ------------------
# comp <- fsl %>%
#   left_join(pf %>% select(fishID, "originalStdLength" = fishStdLength), by = "fishID") %>%
#   rename("recomputedStdLength" = fishStdLength) %>%
#   mutate(diff = recomputedStdLength-originalStdLength,
#          largeDiff = case_when(abs(diff) > 2 ~ T,
#                                TRUE ~ F))
# write.csv(comp, file = here("data", "outputs", "fslComparison.20210320.csv"), row.names = F)
# 
# comp <- comp %>%
#   mutate(diff = abs(recomputedStdLength - originalStdLength))
# 
# comp %>% filter(diff > 2) %>% pull(fishID)
