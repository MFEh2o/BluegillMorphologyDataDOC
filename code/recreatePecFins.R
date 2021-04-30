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
db <- "MFEdb_20210423.db" # name of db

# Load the original pec fin data for comparison ---------------------------
pf <- read.csv(here("data", "unclassified", "PecFinDataNovember_ORIGINAL.csv")) %>%
  select(-imageID)

# Grab FISH_MORPHOMETRICS -------------------------------------------------
fm <- dbTable("fish_morphometrics")

# Grab lakeInfo -----------------------------------------------------------
lakeInfo <- read.csv(here("data", "outputs", "Lake_Info_2020wBasins.csv"))

# Initialize recreated df -------------------------------------------------
pfR <- fm %>%
  select(fishID, imageFile, parameter, parameterValue, replicate) %>%
  filter(replicate == "NA") %>%
  filter(parameter %in% c("pecFinLength", "pecFinBaseWidth", "stdLength")) %>%
  pivot_wider(id_cols = c("fishID", "imageFile"), names_from = "parameter", values_from = "parameterValue")

# Are all the fish represented?
all(pf$fishID %in% pfR$imageFile) # good!

# Limit it to the fish contained in the original file
pfR <- pfR %>%
  filter(imageFile %in% pf$fishID)

nrow(pfR) == nrow(pf)
all(pfR$imageFile == pf$fishID)
all(pfR$pecFinLength == pf$PecFinLengths) # these don't quite match because I recomputed a few of them in the database update
all(pfR$pecFinBaseWidth == pf$baseWidth) # these do match.

# lakeID ------------------------------------------------------------------
pfR <- pfR %>%
  mutate(lakeID = word(fishID, 1, 1, sep = "_")) %>%
  mutate(lakeID = forcats::fct_recode(lakeID,
                                      "Bay" = "BA",
                                      "Birch" = "BH",
                                      "Crampton" = "CR",
                                      "Found" = "FD",
                                      "Hummingbird" = "HB",
                                      "Little_Crooked" = "LC",
                                      "Lost" = "LT",
                                      "McCullough" = "MC",
                                      "Muskellunge" = "MS",
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
              select(lakeName, DOC),
            by = c("lakeID" = "lakeName"))

# Join std lengths from whole-body photos ---------------------------------
l <- fm %>%
  filter(parameter == "stdLength", imageType == "whole_body", replicate == "NA") %>%
  select(fishID, "stdLength" = parameterValue)

nrow(pfR) # 218
pfR <- pfR %>%
  select(-stdLength) %>%
  left_join(l, by = "fishID")
nrow(pfR) # 218

# Size-standardize the pec fins -----------------------------------------
# First, let's compute the mean fishStdLength across all fish
meanFishSize <- mean(pfR$stdLength)

# Pec Fin Lengths
# Now, we make three different models and compare them to determine which coefficient to use. See sizeStandardizationNotes.docx for more details on this approach.
moda <- lm(log(pecFinLength) ~ log(stdLength), data = pfR)
modb <- lm(log(pecFinLength) ~ log(stdLength)+lakeID, data = pfR)
modc <- lm(log(pecFinLength) ~ log(stdLength)*lakeID, data = pfR)

# Compare models to determine which coefficient to use
anova(modc, modb) # model c is significantly better than model b

# We will still use the coefficient from model b, as is standard practice (see sizeStandardizationNotes.docx)
coef <- coef(modb)[2] %>% unname() # pull the "Estimate" parameter from the model
coef # 0.7673796

# Now compute size-standardized fin length, adding the column to `pfR`
## We're using the Kaeuffer version of the formula here (see notes doc)
pfR <- pfR %>%
  mutate(finLengthSS = pecFinLength*(meanFishSize/stdLength)^coef)

## check it against the original
head(pfR$finLengthSS)
head(pf$finLengthSS) # pretty close!

## Pec Fin Base Widths
# Now, we make three different models and compare them to determine which coefficient to use. See sizeStandardizationNotes.docx for more details on this approach.
moda <- lm(log(pecFinBaseWidth) ~ log(stdLength), data = pfR)
modb <- lm(log(pecFinBaseWidth) ~ log(stdLength)+lakeID, data = pfR)
modc <- lm(log(pecFinBaseWidth) ~ log(stdLength)*lakeID, data = pfR)

# Compare models to determine which coefficient to use
anova(modc, modb) # model c is significantly better than model b

# We will still use the coefficient from model b, as is standard practice (see sizeStandardizationNotes.docx)
coef <- coef(modb)[2] %>% unname() # pull the "Estimate" parameter from the model
coef # 0.8631758

# Now compute size-standardized fin base width, adding the column to `pfR`
## We're using the Kaeuffer version of the formula here (see notes doc)
pfR <- pfR %>%
  mutate(finBaseSS = pecFinBaseWidth*(meanFishSize/stdLength)^coef)

## check it against the original
head(pfR$finBaseSS)
head(pf$finBaseSS) # pretty close!

# Compute size-standardized length:width ratio ----------------------------
pfR$finRatioSS <- pfR$finLengthSS/pfR$finBaseSS

## compare with original
head(pfR$finRatioSS)
head(pf$ss_ratio) # similar, but not the same, due to ss body length differences.

# Write out the data ------------------------------------------------------
write.csv(pfR, file = here("data", "outputs", "PecFinDataNovemberFINAL.csv"), row.names = F)
