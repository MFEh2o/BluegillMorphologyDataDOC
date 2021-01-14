#Map
library(maps)
data(lakesMapEnv)
map('usa')
map("state","wisconsin", add=TRUE, fill=TRUE, col='grey', boundary='black')
