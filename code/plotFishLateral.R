# plot_fish_lateral function

# This is taken almost verbatim from the btShapes plotting tutorial, found here: https://aaronolsen.github.io/tutorials/morphometrics/backtransform.html

plot_fish_lateral <- function(xy, coor, size = 1, col = "black"){
  # If 3D, rotate points about x-axis using 3D rotation matrix
  if(ncol(coor) == 3){
    coor <- coor %*% matrix(c(1,0,0,0,cos(-pi/2),sin(-pi/2), 
                              0,-sin(-pi/2),cos(-pi/2)), 
                            nrow=3, ncol=3)
  }
  
  # Get just x,y coordinates (orthographic projection into xy-plane)
  coor <- coor[, 1:2]
  
  # Get plot aspect ratio
  w <- par("pin")[1]/diff(par("usr")[1:2])
  h <- par("pin")[2]/diff(par("usr")[3:4])
  asp <- w/h
  
  # Correct for plot aspect ratio not necessarily being 1:1
  coor[, 1] <- coor[, 1] * (1/asp)
  
  # Scale points and place back in position
  coor <- coor*size
  
  # Center about zero based on range of coordinates
  coor <- coor - matrix(colMeans(apply(coor, 2, range)), 
                        nrow=nrow(coor), ncol=ncol(coor), byrow = TRUE)
  
  # Move shape to PC score
  coor <- coor + matrix(xy, nrow(coor), ncol(coor), byrow = TRUE)
  
  # Set order in which to draw points to create polygon
  polygon_order <- c(1,3:12,16,1) # name landmarks in the order you want them connected
  # Create filled polygon
  polygon(coor[polygon_order, ], col = col, border = col)
}