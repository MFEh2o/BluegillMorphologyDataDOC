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

colFuncLowHigh <- colorRampPalette(c(ltblue, medblue, dkblue, tan, brown))
colFuncHighLow <- colorRampPalette(c(brown, tan, dkblue, medblue, ltblue))
lakeColorsLowHigh <- colFuncLowHigh(14)
lakeColorsHighLow <- colFuncHighLow(14)

# Define shapes for lakes
lakeShapesHighLow <- lakeInfo %>%
  arrange(-DOC) %>%
  mutate(shape = case_when(basin == 4 ~ 21,
                           basin == 7 ~ 22)) %>%
  pull(shape)

lakeShapesLowHigh <- lakeInfo %>%
  arrange(DOC) %>%
  mutate(shape = case_when(basin == 4 ~ 21,
                           basin == 7 ~ 22)) %>%
  pull(shape)

# Function to create a model summary table ---------------------------------
makeSummary <- function(model, conf){
  # Make the table
  table <- data.frame(parameter = c("intercept", "log(DOC)", "basin", 
                                    "log(DOC):basin", "sigmaLake", "sigmaRes", 
                                    "R2marginal", "R2conditional"),
    estimate = c(fixef(model), as.data.frame(VarCorr(model))[1,5], 
                 sigma(model), as.vector(r.squaredGLMM(model))),
    lowerCI = c(conf[c(3:6,1:2),1], NA, NA),
    upperCI = c(conf[c(3:6,1:2),2], NA, NA)
  ) %>%
    mutate(across(c("estimate", "lowerCI", "upperCI"), function(x) round(x, 4))) %>%
    mutate(estimate = case_when(!parameter %in% c("R2marginal", "R2conditional") ~ 
                                  paste0(estimate, " (", lowerCI, ", ", upperCI, ")"),
                                TRUE ~ as.character(estimate))) %>%
    select(-lowerCI, -upperCI)
  
  return(table)
}

# btShapes_wrapper -----------------------------------------------------
# This function plots btShapes with the parameters for the all-individuals PCA. Pulling it out so I don't have to re-type it three times.
btShapes_wrapper <- function(sc, vc, f = plot_fish_lateral, 
                             p = pcs, nn = c(4, 4), mm = dim(landmarks)[2], 
                             rn = dimnames(landmarks)[[1]], 
                             pm = c(0.06, 0.05), sz = 0.038, 
                             cl = gray(0.7)){
  StereoMorph::btShapes(scores = sc, 
                        vectors = vc, 
                        fcn = f, # defined in plotFishLateral.R
                        pcs = p, # the two pc axes for the plot
                        n = nn, # number of shapes to draw along the x and y axes, respectively
                        m = mm, 
                        row.names = rn, # order in which to connect the landmarks
                        pc.margin = pm, 
                        size = sz, col = cl)
}


