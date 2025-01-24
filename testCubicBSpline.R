
setwd("~/COS_LK")

source("bspline_poly.R")
suppressMessages( library( LatticeKrig))

suppressMessages( library( colorspace)) #for some color alternatives

#1d example
dGrid<- seq( -5,5, .001)

look<- CubicBSplineDN(dGrid)
plot( dGrid, look, type="l")

plot( diff(diff( look)), type="l")
plot( diff(diff(diff( look))), type="l")

s<- cbind(rat.diet$t)
z<- cbind(rat.diet$trt)

obj<- LatticeKrig( s, z, NC=10, nlevel=1, a.wght=2.1)
plot( s,z)
lines( s, obj$fitted.values, col="magenta")


LKInfo<- LKrigSetup( s, NC=15, nlevel=1, a.wght=2.1, overlap = 2, NC.buffer = 2, 
                     LKGeometry="LKInterval", BasisFunction = "CubicBSplineDN", 
                     normalize = FALSE)

##new work 1D test
x1 <- cbind(seq(0,105, 0.01))
out <- LKrig.basis(x1, LKinfo = LKInfo)

grid_x <- LKInfo$latticeInfo$grid[[1]]$x

out_full <- spam2full(out)
test_sum <- rowSums(out_full)

plot(x1, test_sum, type = "l", col = "red", ylim = c(0,1.1))
for (j in 1:19) {
  lines(x1, out_full[,j])
}
xline(grid_x)

fit <- LatticeKrig( s, z, LKinfo = LKInfo)

plot( s,z, pch = 16)
lines(s, fit$fitted.values, col = "magenta2")

##new work 2d test
setwd("~/COS_LK")
load(file = "synthetic_data.rda") #sGrid, look_stat synthetic data

synth_df <- data.frame(sGrid, look_stat)

#samples
n <- 500
set.seed(350)
sub1_df <- synth_df[sample(nrow(synth_df), n, replace = FALSE), ]

#sample plot
bubblePlot(sub1_df$x, sub1_df$y, sub1_df$look_stat, 
           col = tim.colors, main = "Synth Data - Sample Points")

#LK setup
s_new <- cbind(sub1_df$x, sub1_df$y)
z_new <- cbind(sub1_df$look_stat)

#cubic bspline
LK_info <- LKrigSetup(s_new, NC=10, nlevel=1, a.wght=4.1, overlap = 2, NC.buffer = 2, 
                      LKGeometry="LKRectangle", BasisFunction = "CubicBSplineDN", 
                      BasisType = "Tensor", normalize = FALSE)

#wendland
LK_infobase <- LKrigSetup(s_new, NC=10, nlevel=1, a.wght=4.1,  
                          LKGeometry="LKRectangle", 
                          BasisType = "Radial", normalize = TRUE)

gridList_new <- list( x= seq(-1, 1, 0.01),
                      y= seq(-1, 1, 0.01))
sGrid_new <- make.surface.grid(gridList_new)

out2d <- LKrig.basis(sGrid_new, LKinfo = LK_info)
out2d_base <- LKrig.basis(sGrid_new, LKinfo = LK_infobase)

test_sum2d <- rowSums(spam2full(out2d))
sum2d_base <- rowSums(spam2full(out2d_base))

#imaging the basis function
set.panel(1,2)
image.plot( as.surface(sGrid_new, test_sum2d), main = "Cubic B-Spline (Tensor)")
rect(-1,-1,1,1, col = NA, border = "magenta", lwd = 2)
image.plot( as.surface(sGrid_new, sum2d_base), main = "Wendland (Radial)")
rect(-1,-1,1,1, col = NA, border = "magenta", lwd = 2)
dev.off()

new_2d <- sum2d_base - test_sum2d
image.plot( as.surface(sGrid_new, new_2d))
#getting LK fits

fit2d <- LatticeKrig( s_new, z_new, LKinfo = LK_info)
fit2d_base <- LatticeKrig( s_new, z_new, LKinfo = LK_infobase)

surface(fit2d, nx = 200, ny = 150, col = tim.colors(256), main = "Cubic B-Spline")
surface(fit2d_base,  nx = 200, ny = 150, col = tim.colors(256), main = "Wendland")

#get difference


M <- 56
gridList_new <- list( x= seq( -1,1,length.out= M),
                      y= seq( -1,1,length.out= M) )

lk_hat <- predictSurface(fit2d, grid.list = gridList_new, nx = 100, ny = 100) #not currently working
image.plot(lk_hat, xlab = "x", ylab = "y", main = "Cubic B-Spline")

lk_hat_base <- predictSurface(fit2d_base, grid.list = gridList_new, nx = 100, ny = 100) #not currently working
image.plot(lk_hat_base, xlab = "x", ylab = "y", main = "Wendland")

lk_diff <- lk_hat
lk_diff$z <- lk_hat$z - lk_hat_base$z
image.plot(lk_diff, xlab = "x", ylab = "y", main = "Diff (B-spline - Wendland)", col = diverge_hsv(256))

set.panel(2,2)
image.plot(lk_hat, xlab = "x", ylab = "y", main = "Cubic B-Spline (Tensor)")
image.plot(lk_hat_base, xlab = "x", ylab = "y", main = "Wendland (Radial)")
bubblePlot(sub1_df$x, sub1_df$y, sub1_df$look_stat, 
           col = tim.colors, main = "Synth Data - Sample Points")
image.plot(lk_diff, xlab = "x", ylab = "y", main = "Diff (B-spline - Wendland)", col = diverge_hsv(256))

#lets explore SE


