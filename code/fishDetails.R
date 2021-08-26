#Additional descriptive stats
#CTS 2021-08-25 for manuscript revisions

# Load packages
library(lubridate)
library(here) # for file paths
source(here("code", "dbUtil.R")) # for the dbTable() function
source(here("code", "defs.R"))
library(RSQLite) # for db connection
library(dplyr)
library(ggplot2)

source(here("code", "defs.R"))

# Connect to the database
dbdir <- here("data") # directory where the database is stored
db <- "MFEdb_20210423.db" # name of db

# Load data
#Load FISH_MORPHOMETRICS
fm <- dbTable("fish_morphometrics")
#Load FISH_INFO
fishInfo <- dbTable("FISH_INFO")
#Convert standardLength to numeric
fishInfo$standardLength <- as.numeric(fishInfo$standardLength)
#Load fishSamples
fishSamples <- dbTable("FISH_SAMPLES")

#Retain fishInfo only for fish that show up in fm
fishInfo <- subset(fishInfo,subset=fishInfo$fishID%in%unique(fm$fishID))

#Check metadataID on these fish
unique(fishInfo$metadataID)  #All 'blgMorphology...' - good.

#Merge in lakeID, dateTimeSet, dateTimeSample, and gear from fishSamples
fishInfo <- merge(fishInfo,fishSamples[,c("sampleID","lakeID","dateTimeSet","dateTimeSample","gear")],by="sampleID",all.x=T,all.y=F)

#Check how many fish we have from each lake 
table(fishInfo$lakeID) #matches table currently in ms - good.

#Summary stats on lengths of fish
lengthMean <- tapply(fishInfo$standardLength,fishInfo$lakeID,mean)
lengthSum <- fishInfo %>%
  group_by(lakeID) %>%
  summarize(lMean=mean(standardLength),
            lMin=min(standardLength),
            lMax=max(standardLength),
            lSD=sd(standardLength))
lengthSum
#Minimum lengths are very low in BA, BH. How many fish that short?
head(sort(fishInfo$standardLength))
#Which fish are these?
shortFish <- fishInfo$fishID[fishInfo$standardLength<15]
shortFish
#When were they sampled?
fishInfo[fishInfo$fishID%in%shortFish,c("fishID","lakeID","dateTimeSample","standardLength")]
#Quick plot of lengths
hist(fishInfo$standardLength)

#Drop the three short fish from fishInfo and fm
fishInfo <- fishInfo[-which(fishInfo$fishID%in%shortFish),]
fm <- fm[-which(fm$fishID%in%shortFish),]

#Re-calculate length summaries
lengthSum <- fishInfo %>%
  group_by(lakeID) %>%
  summarize(lMean=mean(standardLength),
            lMin=min(standardLength),
            lMax=max(standardLength),
            lSD=sd(standardLength))
lengthSum

#Plot length histograms by lake
ggplot(fishInfo) +
  geom_histogram(mapping=aes(standardLength)) +
  facet_wrap(vars(lakeID),ncol=1,strip.position='right') +
  labs(x="Standard length (mm)")