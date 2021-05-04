# Script to recreate eyewidthsFINAL.csv
# Created by Kaija Gahm on 19 February 2021

# Load packages -----------------------------------------------------------
library(dplyr) # for data wrangling
library(tidyr) # for data wrangling
library(stringr) # for text string manipulation
library(ggplot2) # for plotting
library(here) # for file paths
source(here("code", "dbUtil.R")) # for the dbTable() function
source(here("code", "defs.R"))
library(RSQLite) # for db connection

source(here("code", "defs.R"))

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
  filter(replicate == "NA") %>%
  select(lakeID, imageFile, parameter, parameterValue) %>%
  rename("fishID" = imageFile) %>%
  filter(parameter %in% c("eyeWidth", "stdLength")) %>%
  pivot_wider(id_cols = c("lakeID", "fishID"), names_from = "parameter", values_from = "parameterValue")

# Are all the fish represented?
all(ew$fishID %in% ewR$fishID) # good!

# Limit it to the fish contained in the original file
ewR <- ewR %>%
  filter(fishID %in% ew$fishID)

nrow(ewR) == nrow(ew)
all(ewR$fishID == ew$fishID)
all(ewR$eyeWidth == ew$eyeWidth) # all right, these match!

# lakeID ------------------------------------------------------------------
table(ewR$lakeID, exclude = NULL)
table(ew$lakeID, exclude = NULL) # good, the counts line up and there are no NA's. I'm using lake names to avoid incorrect abbreviations.

# Add DOC and basin ------------------------------------------------------
ewR <- ewR %>%
  left_join(lakeInfo %>%
              select(lakeID, DOC, basin),
            by = "lakeID")

sum(is.na(ewR$lakeID)) # no NA's
sum(is.na(ewR$basin)) # no NA's

# Size-standardize the eye widths -----------------------------------------
# First, let's compute the mean fishStdLength across all fish
meanFishSize <- mean(ewR$stdLength)

# Now, we make three different models and compare them to determine which coefficient to use. See sizeStandardizationNotes.docx for more details on this approach.
moda <- lm(log(eyeWidth) ~ log(stdLength), data = ewR)
modb <- lm(log(eyeWidth) ~ log(stdLength)+lakeID, data = ewR)
modc <- lm(log(eyeWidth) ~ log(stdLength)*lakeID, data = ewR)

# Compare models to determine which coefficient to use
anova(modc, modb) # model c is significantly better than model b (p = 0.00229)

# Plot of model C
ewR %>%
  arrange(-DOC) %>% # arrange from highest to lowest DOC
  # reorder the `lakeID` variable from highest to lowest DOC
  mutate(lakeID = factor(lakeID, levels = unique(lakeID))) %>%
  ggplot(aes(x = log(stdLength), y = log(eyeWidth), col = lakeID))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm", alpha = 0.2)+
  scale_color_manual(values = lakeColorsReverse)+
  theme_minimal()+
  labs(title = "Eye width vs body length, log-transformed")
ggsave(filename = "eyeWidthsAllometry.png", path = here("figures", "allometryPlots"), width = 6, height = 4)

# We will still use the coefficient from model b, as is standard practice (see sizeStandardizationNotes.docx)
coef <- coef(modb)[2] %>% unname() # pull the "Estimate" parameter from the model
coef #  0.6266465 

# Now compute size-standardized eye width, adding the column to `ewR`
## We're using the Kaeuffer version of the formula here (see notes doc)
ewR <- ewR %>%
  mutate(eyewidth.ss = eyeWidth*(meanFishSize/stdLength)^coef)

## check it against the original
head(ewR$eyewidth.ss)
head(ew$sizeStandardize) # these are very close but slightly different. I think that's because Chelsea said she rounded the commonWithinGroupSlope parameter to 0.62, whereas I'm using the whole thing. Indeed, I tried it out with the parameter rounded to 0.62, and then they're the same out to like 3 decimals.

# I think that this file is now ready to write out and use in other analyses.

# Write out the data ------------------------------------------------------
write.csv(ewR, file = here("data", "outputs", "eyewidthsFINAL.csv"), row.names = F)
