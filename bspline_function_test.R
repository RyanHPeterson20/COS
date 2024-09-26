setwd("~/COS_LK")
source("DF_LK/R/bspline_poly.R")

#test points
x_grid <- seq( 0.99,6,.01)

min_knot <- floor(min(x_grid))
max_knot <- ceiling(max(x_grid))
buffer <- c(min_knot - 1, max_knot + 1) 

#center of each internal knot
base_knots <- seq(min_knot, max_knot, 1)

#center of each bspline function
x_knots <- seq(min(buffer), max(buffer),1) #including buffer knots


test_bsplines <- matrix(NA, nrow = length(x_grid))
for (j in x_knots) {
  test_bsplines <- cbind(test_bsplines, cubicSplineBasis(x_grid - j))
}
test_bsplines <- test_bsplines[,-1]

#test plot, gray line shows partition of unity
plot(range(x_grid), c(0,1), type = "n", xlab = "x", ylab = "")
abline(v = base_knots, lty = 2)
matlines(x_grid, test_bsplines, ylim = c(0,1), lty = 1)
lines(x_grid, rowSums(test_bsplines), col = "gray", lwd = 2)

#test function as tensor product.

#using the tensor after evaluating
tSpline1 <- outer(look1, look2, "*")
tSpline2 <- outer(look2, look2, "*")

#from LatticeKrig
phi_temp <- Tensor.basis(xyGrid, centered, basis.delta = 1,
                         max.points = NULL,
                         mean.neighbor = 50, 
                         BasisFunction = get(LKinfo$basisInfo$BasisFunction),
                         distance.type = LKinfo$distance.type)

#basis function is derived from LKinfo$basisInfo$BasisFunction:



