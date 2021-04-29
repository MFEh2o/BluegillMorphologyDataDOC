# Script to recreate eyewidthsFINAL.csv
# Created by Kaija Gahm on 19 February 2021

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

# Load the original eyewidthsFINAL for comparison -------------------------
ew <- read.csv(here("data", "unclassified", "eyewidthsFINAL_ORIGINAL.csv")) %>%
  select(-imageID) # this is a file path to chelsea's computer, so don't need this.

# Grab FISH_MORPHOMETRICS -------------------------------------------------
fm <- dbTable("fish_morphometrics")

# Grab lakeInfo -----------------------------------------------------------
lakeInfo <- read.csv(here("data", "outputs", "Lake_Info_2020wBasins.csv"))

# Initialize recreated df -------------------------------------------------
ewR <- fm %>%
  select(imageFile, parameter, parameterValue) %>%
  rename("fishID" = imageFile) %>%
  filter(parameter == "eyeWidth") %>%
  pivot_wider(id_cols = fishID, names_from = "parameter", values_from = "parameterValue")

# Are all the fish represented?
all(ew$fishID %in% ewR$fishID) # good!

# Limit it to the fish contained in the original file
ewR <- ewR %>%
  filter(fishID %in% ew$fishID)

nrow(ewR) == nrow(ew)
all(ewR$fishID == ew$fishID)
all(ewR$eyeWidth == ew$eyeWidth) # all right, these match!

# lakeID ------------------------------------------------------------------
ewR <- ewR %>%
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

table(ewR$lakeID, exclude = NULL)
table(ew$lakeID, exclude = NULL) # good, the counts line up and there are no NA's. I'm using lake names to avoid incorrect abbreviations.

# Add DOC and basin ------------------------------------------------------
ewR <- ewR %>%
  left_join(lakeInfo %>%
              select(lakeName, DOC, basin),
            by = c("lakeID" = "lakeName"))

# Not adding capture method because it doesn't get used in analyses.

# Calculate fish standard length ------------------------------------------
# XXX Soon, will be able to pull this directly from FISH_MORPHOMETRICS
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
  filter(fishID %in% ewR$fishID)

## join
ewR <- ewR %>%
  left_join(fsl, by = "fishID")

## Check
ewR$fishStdLength %>% head(10)
ew$fishStdLength %>% head(10)

# ok, these are the same to 2 decimal places. So not *quite* exactly the same.

# Size-standardize the eye widths -----------------------------------------
# First, let's compute the mean fishStdLength across all fish
meanFishSize <- mean(ewR$fishStdLength)

# Now, we make three different models and compare them to determine which coefficient to use. See sizeStandardizationNotes.docx for more details on this approach.
moda <- lm(log(eyeWidth) ~ log(fishStdLength), data = ewR)
modb <- lm(log(eyeWidth) ~ log(fishStdLength)+lakeID, data = ewR)
modc <- lm(log(eyeWidth) ~ log(fishStdLength)*lakeID, data = ewR)
summary(modc) # we note that there is an interaction effect for CR, LT, and MS.

# Compare models to determine which coefficient to use
anova(modc, modb) # model c is significantly better than model b

# We will still use the coefficient from model b, as is standard practice (see sizeStandardizationNotes.docx)
coef <- coef(modb)[2] %>% unname() # pull the "Estimate" parameter from the model
coef #  0.6266465 

# Now compute size-standardized eye width, adding the column to `ewR`
## We're using the Kaeuffer version of the formula here (see notes doc)
ewR <- ewR %>%
  mutate(eyewidth.ss = eyeWidth*(meanFishSize/fishStdLength)^coef)

## check it against the original
head(ewR$eyewidth.ss)
head(ew$sizeStandardize) # these are very close but slightly different. I think that's because Chelsea said she rounded the commonWithinGroupSlope parameter to 0.62, whereas I'm using the whole thing. Indeed, I tried it out with the parameter rounded to 0.62, and then they're the same out to like 3 decimals.

# I think that this file is now ready to write out and use in other analyses.

# Write out the data ------------------------------------------------------
write.csv(ewR, file = here("data", "outputs", "eyewidthsFINAL.csv"), row.names = F)

