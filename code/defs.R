# Define colors etc. for lakes
brown <- rgb(113, 83, 55, maxColorValue = 255)
tan <- rgb(223, 182, 131, maxColorValue = 255)
dkblue <- rgb(48, 76, 170, maxColorValue = 255)
medblue <- rgb(121, 221, 238, maxColorValue = 255)
ltblue <- rgb(194, 242, 238, maxColorValue = 255)

colfunc <- colorRampPalette(c(ltblue, medblue, dkblue, tan, brown))
lakeColors <- colfunc(14)