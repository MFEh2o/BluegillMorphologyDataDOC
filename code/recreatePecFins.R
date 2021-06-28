# Script to recreate pectoral fin data
# Created by Kaija Gahm on 3 March 2021

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
pfR <- fm %>%
  select(lakeID, fishID, imageFile, parameter, parameterValue, replicate) %>%
  filter(replicate == "NA") %>%
  filter(parameter %in% c("pecFinLength", "pecFinBaseWidth", "stdLength")) %>%
  pivot_wider(id_cols = c("lakeID", "fishID", "imageFile"), 
              names_from = "parameter", values_from = "parameterValue")

# Are all the fish represented?
# Remove any whole-body fishID's (that don't have pec fin measurements)
pfR <- pfR %>%
  filter(!(is.na(pecFinBaseWidth) & is.na(pecFinLength)))

# Add DOC and basin ------------------------------------------------------
pfR <- pfR %>%
  left_join(lakeInfo %>%
              select(lakeID, DOC),
            by = "lakeID")

# Join std lengths from whole-body photos ---------------------------------
l <- fm %>%
  filter(parameter == "stdLength", imageType == "whole_body", replicate == "NA") %>%
  select(fishID, "stdLength" = parameterValue)

nrow(pfR) # 218
pfR <- pfR %>%
  select(-stdLength) %>%
  left_join(l, by = "fishID")
nrow(pfR) # 218 still, good

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

# Plot model C
p <- pfR %>%
  arrange(-DOC) %>% # arrange from highest to lowest DOC
  # reorder the `lakeID` variable from highest to lowest DOC
  mutate(lakeID = factor(lakeID, levels = unique(lakeID))) %>%
  ggplot(aes(x = log(stdLength), y = log(pecFinLength), col = lakeID))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm", alpha = 0.2)+
  scale_color_manual(values = lakeColorsHighLow)+
  theme_minimal()+
  labs(title = "Pec fin length vs body length, log-transformed")
ggsave(p, filename = "pecFinLengthAllometry.png", path = here("figures", "allometryPlots"), width = 6, height = 4)

# We will still use the coefficient from model b, as is standard practice (see sizeStandardizationNotes.docx)
coef <- coef(modb)[2] %>% unname() # pull the "Estimate" parameter from the model
coef # 0.7673796

# Now compute size-standardized fin length, adding the column to `pfR`
## We're using the Kaeuffer version of the formula here (see notes doc)
pfR <- pfR %>%
  mutate(finLengthSS = pecFinLength*(meanFishSize/stdLength)^coef)

## Pec Fin Base Widths
# Now, we make three different models and compare them to determine which coefficient to use. See sizeStandardizationNotes.docx for more details on this approach.
moda <- lm(log(pecFinBaseWidth) ~ log(stdLength), data = pfR)
modb <- lm(log(pecFinBaseWidth) ~ log(stdLength)+lakeID, data = pfR)
modc <- lm(log(pecFinBaseWidth) ~ log(stdLength)*lakeID, data = pfR)

# Compare models to determine which coefficient to use
anova(modc, modb) # model c is significantly better than model b

# Plot model c
p <- pfR %>%
  arrange(-DOC) %>% # arrange from highest to lowest DOC
  # reorder the `lakeID` variable from highest to lowest DOC
  mutate(lakeID = factor(lakeID, levels = unique(lakeID))) %>%
  ggplot(aes(x = log(stdLength), y = log(pecFinBaseWidth), col = lakeID))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm", alpha = 0.2)+
  scale_color_manual(values = lakeColorsHighLow)+
  theme_minimal()+
  labs(title = "Pec fin width vs body length, log-transformed")
ggsave(p, filename = "pecFinWidthAllometry.png", path = here("figures", "allometryPlots"), width = 6, height = 4)

# We will still use the coefficient from model b, as is standard practice (see sizeStandardizationNotes.docx)
coef <- coef(modb)[2] %>% unname() # pull the "Estimate" parameter from the model
coef # 0.8631758

# Now compute size-standardized fin base width, adding the column to `pfR`
## We're using the Kaeuffer version of the formula here (see notes doc)
pfR <- pfR %>%
  mutate(finBaseSS = pecFinBaseWidth*(meanFishSize/stdLength)^coef)

# Compute size-standardized length:width ratio ----------------------------
pfR$finRatioSS <- pfR$finLengthSS/pfR$finBaseSS

# Write out the data ------------------------------------------------------
write.csv(pfR, file = here("data", "outputs", "PecFinDataNovemberFINAL.csv"), row.names = F)