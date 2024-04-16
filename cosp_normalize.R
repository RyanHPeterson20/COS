suppressMessages(library(LatticeKrig))

#set up Lkinfo: (normalize must be false)
M = 200

sDomain<- rbind( c( -1,-1), c( 1,2))
gridDomain <- list(x = seq(-1, 1, length.out = M),
                 y = seq(-1, 2, length.out = M))
domain_grid <- make.surface.grid(gridDomain)

#look up alpha (setting nu as 1 instead)
#a.wght = 4.01 is approx. TPS
LKinfo <- LKrigSetup(sDomain, NC=25, nlevel=1, a.wght = 4.01,
                     NC.buffer = 0, normalize = TRUE, nu = 1)



#load polygroup for testing
#load data 
load("DF_LK/DFSynthData_Base.rda")

#TODO: add in normalize to the basisIntegral.R

#use interp.surface.grid to adjust to the size of the grid for normalize


##temp rebuilding parts of LKrig.basis:
#bulding up some background
nLevel <- LKinfo$nlevel

reducedDomain <- rbind( c( -1,-1), c( 1,1))
gridList <- list(x = seq(-1, 1, length.out = M),
                 y = seq(-1, 1, length.out = M))

boundingBox<- rbind( c( -1,-1),
                     c( 1,-1),
                     c( 1,1),
                     c( -1, 1),
                     c( -1,-1))

xyGrid <- make.surface.grid(gridList)

dx<- gridList$x[2]- gridList$x[1]
dy<- gridList$y[2]- gridList$y[1]

#fix this later with 
zGrid <- Wendland(sqrt(xyGrid[, 1] ^ 2 + xyGrid[, 2] ^ 2),
                  dimension = 2,
                  k = 2)
stdIntegral<- sum( zGrid)*dx*dy


obj <- spam2full( LKrig.basis(domain_grid, LKinfo = LKinfo))
#
N1 <- length(polyGroups)
# total number of basis functions
N2 <- LKinfo$latticeInfo$m
# basis function scales
delta <- LKinfo$latticeInfo$delta

#####left off here (save this viz for latter)
testGrid <- as.data.frame(cbind(xyGrid, zGrid))
bubblePlot(xyGrid, zGrid)
tempGrid <- testGrid[testGrid$zGrid != 0, ]
tempXY <- as.matrix(tempGrid[,c(1,2)])
bubblePlot(tempXY, tempGrid$zGrid)

zlim<- range(zGrid)
surface(as.surface( xyGrid, zGrid) , 
        zlim =zlim, col=turbo(256))

#directly from LKrig.basis

l = 1
basis.delta <- LKrigLatticeScales(LKinfo) 
centers<- LKrigLatticeCenters( LKinfo, Level=l ) #add for l in 1:nlevel loop
PHItemp <- Radial.basis(  domain_grid, centers, basis.delta[l],
                          max.points = LKinfo$basisInfo$max.points,
                          mean.neighbor = LKinfo$basisInfo$mean.neighbor, 
                          BasisFunction = get(LKinfo$basisInfo$BasisFunction),
                          distance.type = LKinfo$distance.type,
                          verbose = FALSE)

#TODO:  LKrigSAR.LKRectangle if(first.order) has condition with length > 1
##- figure out why there is an issue
#LKrigSAR(LKinfo, Level = Level)
#first.order<- attr( object$a.wght, "first.order")[Level]
attr( LKinfo$a.wght, "first.order")[1]

#change below to test_wght
test_wght <- LKrigNormalizeBasis(LKinfo,  Level=l, PHItemp)   #working at the moment
  
#above issue with LKrigSAR.LKRectangle if(first.order) has condition with length > 1
##- figure out why there is an issue
#S3method(LKrigSAR,LKRectangle)

test_normalize <- 1/sqrt(test_wght)

test_wght_fast <- LKrigNormalizeBasisFast(LKinfo,  Level=1,  x=xyGrid)   

#viz for normalize wghts
bubblePlot(domain_grid, test_normalize)



#adjust test_normalize down to the the constrained domain
wght_list <- list(x = domain_grid[,1], y = domain_grid[,2], z= test_normalize)
xy_list <- list(x = xyGrid[,1], y = xyGrid[,2])
bounding_normalize <- interp.surface(wght_list, xyGrid)

#interp.surface.grid code block
obj <- wght_list
grid.list <-xy_list
x <- grid.list$x
y <- grid.list$y
M <- length(x)
N <- length(y)
out <- matrix(NA, nrow = M, ncol = N)
for (i in 1:M) {
  out[i, ] <- interp.surface(obj, cbind(rep(x[i], N), y))
}
list(x = x, y = y, z = out)

i <- 1
test_x <- interp.surface(obj, cbind(rep(x[i], N), y))
#now we have test_normalize and zGrid with the same length


#TODO: create viz as we apply the basis integration


#TODO: pick up where I left off, and double check how to multiply the wghts--