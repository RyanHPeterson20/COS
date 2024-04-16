
suppressMessages(library(LatticeKrig))

#for adding in choice of basis types in the basisIntegral function

#GOAL:
##1. Show that we can produce a radial.basis that is equivalent to zGrid
##2. Repeat the process with tensor.basis for basisType options

#TODO: add everything up to the the zGrid
M <- 100
gridList <- list(x = seq(-1, 1, length.out = M),
                 y = seq(-1, 1, length.out = M))

xyGrid <- make.surface.grid(gridList)

zGrid <- Wendland(sqrt(xyGrid[, 1] ^ 2 + xyGrid[, 2] ^ 2),
                  dimension = 2,
                  k = 2)

##including LKinfo
#full domain
sDomain<- rbind( c( -2,-2), c( 2,2))
gridDomain <- list(x = seq(-2, 2, length.out = M),
                   y = seq(-2, 2, length.out = M))
domain_grid <- make.surface.grid(gridDomain) #not needed at the moment

#a.wght = 4.01 is approx. TPS, nu = 1
LKinfo <- LKrigSetup(sDomain, NC=5, nlevel=1, a.wght = 4.01,
                     NC.buffer = 0, normalize = FALSE, nu = 1)



#testing Radial.basis

LKinfo_radial <- LKrigSetup(sDomain, NC=5, nlevel=1, a.wght = 4.01,
                            NC.buffer = 0, normalize = FALSE, nu = 1,
                            BasisType = "Radial")

#note: mean.neighbors can be 1
centered <- matrix(c(0,0), nrow = 1)
colnames(centered) <- c("x", "y")
phi_test <- Radial.basis(xyGrid, centered, basis.delta = 1,
                         max.points = NULL,
                         mean.neighbor = 50, 
                         BasisFunction = get(LKinfo_radial$basisInfo$BasisFunction),
                         distance.type = LKinfo_radial$distance.type,
                         verbose = FALSE)

test_phi <- spam2full(phi_test)
test.for.zero(zGrid, test_phi)

#testing: tensor.basis



#------------------test section-------------------------#

#Full radial basis method from the lkinfo 
l = 1 #add for l in 1:nlevel loop
basis.delta <- LKrigLatticeScales(LKinfo) 
centers <- LKrigLatticeCenters( LKinfo, Level=1 ) 
PHItemp <- Radial.basis(  domain_grid, centers, basis.delta[l],
                          max.points = LKinfo$basisInfo$max.points,
                          mean.neighbor = LKinfo$basisInfo$mean.neighbor, 
                          BasisFunction = get(LKinfo$basisInfo$BasisFunction),
                          distance.type = LKinfo$distance.type,
                          verbose = FALSE)

#test lk distance function
out <- LKrigDistance(xyGrid, centered, basis.delta = 1,
                      max.points = NULL,
                      mean.neighbor = 50, 
                      distance.type = distance.type,
                      components = FALSE)
x1 <- xyGrid
x2 <- centered
n1 <- nrow(x1)
n2 <- nrow(x2)
dimension<- ncol(x1)
Nmax <- max(c(n1,n2)) * 50
