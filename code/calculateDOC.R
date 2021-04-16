# Script to re-create DOC values by pulling data from the database
# Created by Kaija Gahm on 15 April 2021

# Load packages -----------------------------------------------------------
library(tidyverse) # for data wrangling
library(here) # for file paths
library(RSQLite) # for database connection

# Connect to the database -------------------------------------------------
con <- dbConnect(SQLite(), here("data", "MFEdb_20210402.db"))

# List of lakes -----------------------------------------------------------
l <- read.csv(here("data", "inputs", "lakeInfo.csv"))
lakes <- l$lakeID

# Pull in the DOC data from WATER_CHEM ------------------------------------
wc <- tbl(con, "water_chem") %>%
  collect()

wc2 <- wc %>%
  filter(parameter == "DOC",
         lakeID %in% lakes) %>%
  collect()

# Get only PML data,  and summarize by lake -------------------------------
pml_doc <- wc2 %>%
  mutate(parameterValue = as.numeric(parameterValue)/1000) %>% # divide by 1000 to convert from ug/L to mg/L.
  filter(grepl("PML", sampleID), # only PML data
         flag == "0", # remove any flagged data
         dateSample < "2019-01-01") %>% # take all data up through 2018
  group_by(lakeID) %>%
  summarize(meanDOC = mean(parameterValue, na.rm = T),
            sdDOC = sd(parameterValue, na.rm = T),
            minDOC = min(parameterValue, na.rm = T),
            maxDOC = max(parameterValue, na.rm = T),
            medDOC = median(parameterValue, na.rm = T))

# Join calculated values with manually-selected values --------------------
l2 <- l %>%
  select(lakeID, DOC) %>%
  left_join(pml_doc, by = "lakeID") # oh, that's weird--not all of these joined properly.

# Do we maybe not have samples for the other lakes?
wc %>% filter(lakeID == "FD", 
              grepl("PML", sampleID), 
              parameter == "DOC") # oh that's weird...

wc %>% filter(lakeID == "FD",
              parameter == "DOC") # we only have one DOC datum from this lake?

## Let's get a full summary of which data we have, and from which lakes
### First, check that all the lakeID's are correct.
lakes
lakesdb <- dbReadTable(con, "lakes")
all(lakes %in% lakesdb$lakeID) # okay, good, they're all in here. And I checked manually that they're all the right ones.

years <- 2010:2018
docData <- lapply(years, function(x) wc %>% 
                    filter(lakeID %in% lakes) %>%
                    filter(lubridate::year(dateSample) == x) %>%
                    filter(parameter == "DOC") %>%
                    mutate(depthClass = word(sampleID, 5, 5, "_")) %>%
                    group_by(lakeID, depthClass) %>%
                    summarize(n = n()) %>%
                    pivot_wider(id_cols = lakeID, names_from = depthClass, values_from = n))
docData

## Looks like all lakes had at least one surface sample in 2018. Let's look at those.
wc %>%
  filter(lakeID %in% lakes, parameter == "DOC") %>%
  mutate(depthClass = word(sampleID, 5, 5, "_")) %>%
  filter(depthClass == "surface",
         lubridate::year(dateSample) == "2018") %>%
  group_by(lakeID) %>%
  summarize(minDate = min(dateSample),
            maxDate = max(dateSample),
            n = n()) # all have aug/september samples

## The closest thing we have for BA is an august 2015 PML sample. Chris says that's close enough.
wc %>% filter(lakeID == "BA", parameter == "DOC", grepl("-08-|-09-", dateSample)) %>% View()

final <- wc %>%
  filter(lakeID %in% lakes, parameter == "DOC") %>%
  mutate(depthClass = word(sampleID, 5, 5, "_")) %>%
  filter(depthClass == "surface",
         lubridate::year(dateSample) == "2018",
         lubridate::month(dateSample) %in% c("8", "9")) %>%
  arrange(lakeID, dateSample) %>%
  group_by(lakeID) %>%
  slice(1) %>%
  ungroup()# take the earlier sample, for the lakes that have more than one (i.e. early august, as opposed to late september)

baFinal <- wc %>% filter(lakeID == "BA", parameter == "DOC", grepl("-08-|-09-", dateSample), grepl("PML", sampleID), lubridate::year(dateSample) == "2015")

final <- final %>%
  bind_rows(baFinal)

# Plot against original values --------------------------------------------
## Bind on the values from lakeInfo
final <- final %>%
  mutate(parameterValue = as.numeric(parameterValue)/1000) %>%
  left_join(l %>% select(lakeID, DOC), by = "lakeID") %>%
  rename(docOrig = DOC,
         docDB = parameterValue)

## Plot
final %>%
  ggplot(aes(x = docOrig, y = docDB))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_abline(slope = 1)

# Instead,  a plot averaging all DOC values from surface and PML --------
all <- wc %>%
  filter(lakeID %in% lakes, parameter == "DOC", grepl("PML|surface", sampleID), lubridate::year(dateSample) < 2019) %>%
  mutate(parameterValue = as.numeric(parameterValue)/1000) %>%
  filter(flag == "0") %>%
  left_join(l %>% select(lakeID, DOC), by = "lakeID")


## plot all data
all %>%
  mutate(lakeID = fct_reorder(.f = lakeID, .x = parameterValue, .desc = T)) %>%
  ggplot(aes(x = lakeID, y = parameterValue))+
  geom_jitter(aes(col = lubridate::year(dateSample)), width = 0.2, alpha = 0.8)+
  scale_color_viridis_c(begin = 0, end = 0.9)+
  geom_point(aes(x = lakeID, y = DOC), col = "red")

## plot averages
all %>%
  group_by(lakeID) %>%
  summarize(docDB_avgall = mean(parameterValue)) %>%
  left_join(l %>% select(lakeID, DOC), by = "lakeID") %>%
  rename("docOrig" = DOC) %>%
  ggplot(aes(x = docOrig, y = docDB_avgall))+
  geom_point()+
  geom_smooth(method = "lm")+
  geom_abline(slope = 1)



