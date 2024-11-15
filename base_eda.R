
#library
suppressMessages(library(LatticeKrig))
suppressMessages(library(fields))
suppressMessages(library(spam64))

suppressMessages(library(tictoc))
suppressMessages(library(MCMCpack))


#load in data via:

#TODO: double check the data, and convert the polygons to centroid then plot (via bubbleplot)
setwd("~/COS_LK/Bayes")

#TODO: figure out which .rda is the correct data to use, currently testing these
load("Data/trop_data_full.rda")
load("Data/mop_data_full.rda")


#get centroids for each, then constrain to SE Aus for simplicity.

##centroid function
polygon_centroid <- function(polygon) {
  # Ensure the polygon closes by repeating the first point at the end
  if (!all(polygon[1, ] == polygon[nrow(polygon), ])) {
    polygon <- rbind(polygon, polygon[1, ])
  }
  
  # Extract x and y coordinates
  x <- polygon[, 1]
  y <- polygon[, 2]
  
  # Calculate the area A
  A <- 0.5 * sum(x[-1] * y[-length(y)] - x[-length(x)] * y[-1])
  
  # Calculate C_x and C_y
  C_x <- (1 / (6 * A)) * sum((x[-1] + x[-length(x)]) * (x[-1] * y[-length(y)] - x[-length(x)] * y[-1]))
  C_y <- (1 / (6 * A)) * sum((y[-1] + y[-length(y)]) * (x[-1] * y[-length(y)] - x[-length(x)] * y[-1]))
  
  # Return centroid as a vector (change this to 2 column matrix)
  return(matrix(c(x = C_x, y = C_y), ncol = 2))
}

## bounding with Eastern Aus (NE/SE Aus)
#NE aus 134 : 155 long, -25 : -10 lat
#SE aus 134 : 155 long, -48 : -25 lat

##example: list(c(134,-10), c(155, -10), c(155, -48), c(134, -48), c(134,-10))

aus_bound <- rbind(c(134, -10), 
                   c(155, -10), 
                   c(155, -48), 
                   c(134, -48), 
                   c(134, -10))
colnames(aus_bound) <- c("x", "y")

NEaus_bound <- rbind(c(134, -10), 
                     c(155, -10), 
                     c(155, -25), 
                     c(134, -25), 
                     c(134, -10))
colnames(NEaus_bound) <- c("x", "y")

SEaus_bound <- rbind(c(134, -25), 
                     c(155, -25), 
                     c(155, -48), 
                     c(134, -48), 
                     c(134, -25))
colnames(SEaus_bound) <- c("x", "y")


aus_polycentroid <- matrix(NA, ncol = 3)
colnames(aus_polycentroid) <- c("x", "y", "z")

NEaus_polycentroid <- matrix(NA, ncol = 3)
colnames(NEaus_polycentroid) <- c("x", "y", "z")

SEaus_polycentroid <- matrix(NA, ncol = 3)
colnames(SEaus_polycentroid) <- c("x", "y", "z")

for (i in 1:length(polyGroups_mopFull)) {
  
  temp_point <- polygon_centroid(polyGroups_mopFull[[i]])
  temp_data <- mop1_coFull[i]
  
  if (fields::in.poly( temp_point, aus_bound)) {
      aus_polycentroid <- rbind(aus_polycentroid, c(temp_point, temp_data ))
      
      if (fields::in.poly( temp_point, NEaus_bound)) {
        NEaus_polycentroid <- rbind(NEaus_polycentroid, c(temp_point, temp_data ))
      }
      
      if (fields::in.poly( temp_point, SEaus_bound)) {
        SEaus_polycentroid <- rbind(SEaus_polycentroid, c(temp_point, temp_data ))
      }
      
  }
  
}

aus_polycentroid <- aus_polycentroid[-1, ]
NEaus_polycentroid <- NEaus_polycentroid[-1, ]
SEaus_polycentroid <- SEaus_polycentroid[-1, ]

aus_df <- as.data.frame(aus_polycentroid)
NEaus_df <- as.data.frame(NEaus_polycentroid)
SEaus_df <- as.data.frame(SEaus_polycentroid)

#adding in log response
NEaus_df$log_z <- log(NEaus_df$z)
SEaus_df$log_z <- log(SEaus_df$z)



bubblePlot(aus_df$x, aus_df$y, log(aus_df$z), col = tim.colors)
world(add = TRUE)

bubblePlot(NEaus_df$x, NEaus_df$y, log(NEaus_df$z), col = tim.colors)
world(add = TRUE)

bubblePlot(SEaus_df$x, SEaus_df$y, log(SEaus_df$z), col = tim.colors)
world(add = TRUE)

png(filename = "mop_neaus.png", width = 1800, height = 1800, res = 300)
bubblePlot(NEaus_df$x, NEaus_df$y, log(NEaus_df$z), col = tim.colors, main = "MOPITT - Total Column CO", 
           xlab = "Longitude", ylab = "Latitude")
mtext("log(mol/cm^2)", side = 1, adj = 1.4, cex = 1)
world(add = TRUE, lwd = 2)
dev.off()

png(filename = "mop_seaus.png", width = 1800, height = 1800, res = 300)
bubblePlot(SEaus_df$x, SEaus_df$y, log(SEaus_df$z), col = tim.colors, main = "MOPITT - Total Column CO", 
           xlab = "Longitude", ylab = "Latitude")
mtext("log(mol/cm^2)", side = 1, adj = 1.4, cex = 1)
world(add = TRUE, lwd = 2)
dev.off()

## check vgrams
SE_coords <- cbind(SEaus_df$x, SEaus_df$y)
variogram_SEaus <- vgram(SE_coords, SEaus_df$log_z, N = 20, dmax = 15)

plot(variogram_SEaus, main = "SE Aus - Empirical Variogram", xlab = "Distance", ylab = "Semivariance")

NE_coords <- cbind(NEaus_df$x, NEaus_df$y)
variogram_NEaus <- vgram(NE_coords, NEaus_df$log_z,  N = 20, dmax = 15)

plot(variogram_NEaus, main = "NE Aus - Empirical Variogram", xlab = "Distance", ylab = "Semivariance")



#mop_centroid <- polygon_centroid(polyGroups_mopFull)
##fields::in.poly(  polyStd, boundingBox)

# fields MLE
#NOTE: these take forever
tic()
fit_SE <- spatialProcess(SE_coords, SEaus_df$log_z)
toc()

tic()
fit_NE <- spatialProcess(NE_coords, NEaus_df$log_z)
toc()

#TODO: add 
fhat_SE <- predictSurface(fit_SE)
image.plot(fhat_SE)
world(add = TRUE)

fhat_NE <- predictSurface(fit_NE)
image.plot(fhat_NE)
world(add = TRUE)


summary(fit_SE)


#artificial data

library(LatticeKrig)
library(spam64)

# DISCLAIMER: This is just scratching the surface. 
# To produce spatially varying sigma2, spatially varying awght
# or even anisotropic awght, check out the LK examples here: 

# Setting up the spatial domain
sDomain <- cbind(c(-1,1),
                 c(-1,1))
M <- 96
gridList<- list( x= seq( -1,1,length.out= M),
                 y= seq( -1,1,length.out= M) )
sGrid <- make.surface.grid(gridList)

######### How to simulate with single awght ######### 
LKinfo_stat <- LKrigSetup(sDomain, nu=1.0, nlevel=1, NC.buffer = 2,
                          a.wght = 4.01,  lambda = 0.5, sigma = 1, rho = 2,
                          NC=20, normalize=FALSE)

set.seed(351)
look_stat <- LKrig.sim(sGrid, LKinfo_stat, M=1)

#TODO: combine sGrid and look_stat at some point
save(sGrid, look_stat, LKinfo_stat, file = "synthetic_data.rda")

png(filename = "synth.png", width = 1800, height = 1800, res = 300)
image.plot(matrix(look_stat,96,96), 
           xlab = "x", ylab = "y")
dev.off()

shape <- 3   # shape parameter α
scale <- 1   # scale parameter β

# Plot the inverse-gamma density
curve(dinvgamma(x, shape = shape, scale = scale),
      from = 0.1, to = 5,   # Define the x-axis range, avoiding zero
      ylab = "Density",
      xlab = "x",
      xlim = c(0,5),
      main = paste("Inverse-Gamma Distribution\n(shape =", shape, ", scale =", scale, ")"),
      col = "blue",
      lwd = 2)


