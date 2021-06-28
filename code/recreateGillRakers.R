# Script to recreate Gill_Rakers_2018_Final.csv
# Created by Kaija Gahm on 25 February 2021

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
grR <- fm %>%
  select(lakeID, replicate, imageFile, parameter, parameterValue) %>%
  rename("fishID" = imageFile)

## filter to include only the gill raker parameters, and recode them to match the column names in the original Gill_Rakers_2018 sheet
grR <- grR %>%
  filter(replicate == "NA") %>%
  filter(grepl("Raker", parameter)|parameter == "stdLength") %>%
  mutate(parameter = str_replace(parameter, "length", "length_"),
         parameter = str_replace(parameter, "spaceRaker", "space")) %>%
  pivot_wider(id_cols = c("fishID", "lakeID"),
              names_from = "parameter",
              values_from = "parameterValue") %>%
  rename("total_RakerNum" = "totalRakerCount")

# Remove any fish that don't have raker measurements
grR <- grR %>%
  filter(!is.na(length_Raker1))

# Add DOC and basin ------------------------------------------------------
grR <- grR %>%
  left_join(lakeInfo %>%
              select(lakeID, "lakeDOC" = DOC, basin),
            by = "lakeID")
sum(is.na(grR$lakeID))
sum(is.na(grR$basin))

# Compute average raker lengths and spaces --------------------------------
grR <- grR %>%
  rowwise() %>%
  mutate(avgRakerLength = mean(c_across(contains("length_Raker"))),
         avgRakerSpace = mean(c_across(contains("space"))),
         avgL_4.7 = mean(c_across(c(length_Raker4, length_Raker5, length_Raker6, length_Raker7))),
         avgS_4.6 = mean(c_across(c(space4, space5, space6))))

# Size-standardize the averages -------------------------------------------
# The only size-standardized measurements that are used in the final analysis are the avgL_4.7 measurement and the avgS_4.6 measurement, both size-standardized independent of lakeID (so, size-standardized using the full data set, as seen in recreateEyewidths.R).
# So, I'm only going to perform size-standardizations for those columns.

# Average Raker Length (Rakers 4-7)
# First, let's compute the mean fishStdLength across all fish
meanFishSize <- mean(grR$stdLength) 

# Now, we make three different models and compare them to determine which coefficient to use. See sizeStandardizationNotes.docx for more details on this approach.
moda <- lm(log(avgL_4.7) ~ log(stdLength), data = grR)
modb <- lm(log(avgL_4.7) ~ log(stdLength)+lakeID, data = grR)
modc <- lm(log(avgL_4.7) ~ log(stdLength)*lakeID, data = grR)

# Compare models to determine which coefficient to use
anova(modc, modb) # model c is significantly better than model b

# Plot of model C
p <- grR %>%
  arrange(-lakeDOC) %>% # arrange from highest to lowest DOC
  # reorder the `lakeID` variable from highest to lowest DOC
  mutate(lakeID = factor(lakeID, levels = unique(lakeID))) %>%
  ggplot(aes(x = log(stdLength), y = log(avgL_4.7), col = lakeID))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm", alpha = 0.2)+
  scale_color_manual(values = lakeColorsHighLow)+
  theme_minimal()+
  labs(title = "Raker length vs body length, log-transformed")

ggsave(p, filename = "rakerLengthAllometry.png",
       path = here("figures", "allometryPlots"), 
       width = 6, height = 4)

# We will still use the coefficient from model b, as is standard practice (see sizeStandardizationNotes.docx)
coef <- coef(modb)[2] %>% unname() # pull the "Estimate" parameter from the model
coef # 0.398793

# Now compute size-standardized raker length, adding the column to `grR`
## We're using the Kaeuffer version of the formula here (see notes doc)
grR <- grR %>%
  mutate(avgL2_ss = avgL_4.7*(meanFishSize/stdLength)^coef)

# Same thing for the Average Raker Space (Spaces 4-6)
# Make three different models and compare them to determine which coefficient to use. See sizeStandardizationNotes.docx for more details on this approach.
moda <- lm(log(avgS_4.6) ~ log(stdLength), data = grR)
modb <- lm(log(avgS_4.6) ~ log(stdLength)+lakeID, data = grR)
modc <- lm(log(avgS_4.6) ~ log(stdLength)*lakeID, data = grR)

# Compare models to determine which coefficient to use
anova(modc, modb) # model c is significantly better than model b

# Plot of model C
p <- grR %>%
  arrange(-lakeDOC) %>% # arrange from highest to lowest DOC
  # reorder the `lakeID` variable from highest to lowest DOC
  mutate(lakeID = factor(lakeID, levels = unique(lakeID))) %>%
  ggplot(aes(x = log(stdLength), y = log(avgS_4.6), col = lakeID))+
  geom_point(alpha = 0.5)+
  geom_smooth(method = "lm", alpha = 0.2)+
  scale_color_manual(values = lakeColorsHighLow)+
  theme_minimal()+
  labs(title = "Raker space vs body length, log-transformed")
ggsave(p, filename = "rakerSpaceAllometry.png", 
       path = here("figures", "allometryPlots"), width = 6, height = 4)

# We will still use the coefficient from model b, as is standard practice (see sizeStandardizationNotes.docx)
coef <- coef(modb)[2] %>% unname() # pull the "Estimate" parameter from the model
coef # 0.72846

# Now compute size-standardized raker space, adding the column to `grR`
## We're using the Kaeuffer version of the formula here (see notes doc)
grR <- grR %>%
  mutate(avgS2_ss = avgS_4.6*(meanFishSize/stdLength)^coef)

# Write out the data ------------------------------------------------------
write.csv(grR, 
          file = here("data", "outputs", "Gill_Rakers_2018_Final.csv"), 
          row.names = F)