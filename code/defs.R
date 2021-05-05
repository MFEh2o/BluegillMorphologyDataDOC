# Load packages
library(grDevices)
library(here)

# Read in the lake info so we can assign colors by DOC value
lakeInfo <- read.csv(here("data", "outputs", "Lake_Info_2020wBasins.csv"))

# Get lakes ordered by DOC, lowest to highest
lakesLowHigh <- lakeInfo %>%
  arrange(DOC) %>%
  pull(lakeID)

# Get lakes ordered by DOC, highest to lowest
lakesHighLow <- lakeInfo %>%
  arrange(-DOC) %>%
  pull(lakeID)

# Define colors etc. for lakes
brown <- rgb(113, 83, 55, maxColorValue = 255)
tan <- rgb(223, 182, 131, maxColorValue = 255)
dkblue <- rgb(48, 76, 170, maxColorValue = 255)
medblue <- rgb(121, 221, 238, maxColorValue = 255)
ltblue <- rgb(194, 242, 238, maxColorValue = 255)

# XXX rename these to specify low to high or high to low, and change the rest of the code accordingly.
colfunc <- colorRampPalette(c(ltblue, medblue, dkblue, tan, brown))
colfuncreverse <- colorRampPalette(c(brown, tan, dkblue, medblue, ltblue))
lakeColors <- colfunc(14)
lakeColorsReverse <- colfuncreverse(14)
