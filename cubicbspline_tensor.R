suppressMessages(library(LatticeKrig))

#setwd("~/COS_LK")

source("DF_LK/R/bspline_poly.R")


#2d example
x_grid <- seq(0, 6, 0.01)
y_grid <- seq(0, 6, 0.01)
nc <- 8 #number of grid points
buffer <- 1
overlap <- 1  #fixed

range_xgrid <-  range(x_grid)
range_ygrid <-  range(y_grid)

grid.info_x <- list( xmin = range_xgrid[1], xmax= range_xgrid[2], 
                     range = rbind(range_xgrid))
grid.info_y <- list( ymin = range_ygrid[1], ymax= range_ygrid[2], 
                     range = rbind(range_ygrid))

delta <- ( grid.info_x$xmax - grid.info_x$xmin ) / ( nc - 1 )
buffer.width <- buffer * delta

#lattice centers
grid.list_full <- list(x = seq(grid.info_x$xmin - buffer.width, 
                               grid.info_x$xmax + buffer.width, delta),
                       y = seq(grid.info_y$ymin - buffer.width, 
                               grid.info_y$ymax + buffer.width, delta)) #centers

class( grid.list_full) <- "gridList"  

basis.delta <- (delta*overlap)

#b-spline test
test.grid_x <- seq(min(grid.list_full$x-1), max(grid.list_full$x+1), 0.01)
test.grid_y <- seq(min(grid.list_full$y-1), max(grid.list_full$y+1), 0.01)

test_bsplines_x <- matrix(NA, nrow = length(test.grid_x))
for(j in grid.list_full$x){
  test_bsplines_x <- cbind(test_bsplines_x, cubicSplineBasis(((test.grid_x-j)/basis.delta)))
}

test_bsplines_y <- matrix(NA, nrow = length(test.grid_y))
for(k in grid.list_full$y){
  test_bsplines_y <- cbind(test_bsplines_y, cubicSplineBasis(((test.grid_y-k)/basis.delta)))
}

test_bsplines_x <- test_bsplines_x[ ,-1]
test_bsplines_y <- test_bsplines_y[ ,-1]

tensor_list <- list()
for (i in 1:length(grid.list_full$x)) {
  for (k in 1:length(grid.list_full$y)) {
    phi_x <- test_bsplines_x[ ,i]
    phi_y <- test_bsplines_y[ ,k]
    tensor_xy <- outer(phi_x, phi_y, "*")
    tensor_list[[paste0("x",i,"y",k)]] <- tensor_xy
  }
}

tensor_sum <-  Reduce("+", tensor_list)

xy_surface <- make.surface.grid( list(x=test.grid_x, y= test.grid_y) )
image.plot( as.surface(xy_surface, tensor_sum))
rect(0,0,6,6, col = NA, border = "magenta", lwd = 2)
