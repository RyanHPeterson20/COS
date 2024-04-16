#looking into each element part of the basisIntegral function.

#TODO:
## -define choice for alpha or nu when using nLevel > 1
## -finish up replicating the Basis integral function (basisIntegral)
## -adjust the bounding box to accurately represent the compact support of the function 
### --Check the standardized polygons first to make sure we need to change this from box to radial.

#parameters
##polyGroups,
##LKinfo,
##M = 200, 
##check=TRUE

suppressMessages(library(LatticeKrig))

#set up Lkinfo: (normalize must be false)
sDomain<- rbind( c( -1,-1), c( 1,2))

#look up alpha (setting nu as 1 instead)
LKinfo <- LKrigSetup(sDomain, NC=25, nlevel=2, a.wght = 4.1,
                    NC.buffer = 0, normalize = FALSE, nu = 1)
#check alpha setup for multiple levels, missing alpha specification
##see paper for more details on alpha selection.

M = 200

#load polygroup for testing
#load data 
load("DF_LK/DFSynthData_Base.rda")


#-------------------------------------------------#
#isolated testing of function, run each part seperately:
polyGroups <- polyGroups #try polyGroupsWide for actual test


nLevel <- LKinfo$nlevel



gridList <- list(x = seq(-1, 1, length.out = M),
                 y = seq(-1, 1, length.out = M))

#what does this do? Can we change this to a round object
boundingBox<- rbind( c( -1,-1),
                     c( 1,-1),
                     c( 1,1),
                     c( -1, 1),
                     c( -1,-1))

dx<- gridList$x[2]- gridList$x[1]
dy<- gridList$y[2]- gridList$y[1]


xyGrid <- make.surface.grid(gridList)

# values for the standard basis function on the grid
#Note: the distance is defined as the L2 norm.
#TODO: change this reflect the covariance set up in LKinfo
#use LKrig.basis for basis functions
test_basis <- LKrig.basis(xyGrid, LKinfo, verbose = TRUE) #in spam (get back from spam or use as is)
test_basis <- spam2full(test_basis)
zGrid <- Wendland(sqrt(xyGrid[, 1] ^ 2 + xyGrid[, 2] ^ 2),
                  dimension = 2,
                  k = 2)

#test plot of zGrid 
desc_zGrid <- sort(zGrid, decreasing = TRUE)
testGrid <- as.data.frame(cbind(xyGrid, zGrid))
bubblePlot(xyGrid, zGrid)


#select for values not equal to zero
tempGrid <- testGrid[testGrid$zGrid != 0, ]
tempXY <- as.matrix(tempGrid[,c(1,2)])
bubblePlot(tempXY, tempGrid$zGrid)

#TODO: create a smooth surface to plot this


#just an integral over a wendland basis function of d \in [0,1]
stdIntegral<- sum( zGrid)*dx*dy


#------ Integral set-up---##
#not verified yet, copy to basisIntegralNew.R when done

# total number of data "blocks"
N1 <- length(polyGroups)
# total number of basis functions 
N2 <- LKinfo$latticeInfo$m
#TODO: make adjusments for level
test_level <- LKinfo$latticeInfo$mLevel
#potential change:
##N2 <- sum(LKinfo$latticeInfo$mLevel)

# basis function scales (distance between lattice points?)
delta <- LKinfo$latticeInfo$delta
integral <- NULL
J <- NULL
K <- NULL

checkIntegral<- NULL # stays NULL if check == FALSE

#parallelize here
#add in libraries for parallel:
#TODO: figure out what libraries are absolutely necesssary
suppressMessages(library(foreach))
suppressMessages(library(parallel))
suppressMessages(library(doParallel))

cores <- 8

#cluster set-up
n.cores <- cores

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

#TODO: fix the the combine method and correct output
basisInt <- foreach(j = 1:N1, .combine=rbind) %dopar% {
  polyTmp <- polyGroups[[j]]
  
  for (l in 1:nLevel) {
    # get info for the L^th level of the multi-resolution basis. 
    basisIndex <-(1:LKinfo$latticeInfo$mLevel[l])
    basisOffset <- LKinfo$latticeInfo$offset[l]
    basisScale <- LKinfo$latticeInfo$delta[l]*LKinfo$basisInfo$overlap
    basisGridList <- LKinfo$latticeInfo$grid[[l]]
    basisCenters <- fields::make.surface.grid(basisGridList) #requires
    
    tempIntegral <- NULL
    for (k in basisIndex){
      # k is position at the L^th level
      polyStd <-
        cbind((polyTmp[, 1] - basisCenters[k, 1]) / basisScale,
              (polyTmp[, 2] - basisCenters[k, 2]) / basisScale)
      
      allOutside<- fields::in.poly(  polyStd, boundingBox)
      # only look at polygonal if at least one point is inside the bounding box
      # note that integrals with basis function support outside of the
      # polygon are skipped and by defualt set to zero
      
      if(sum(allOutside) > 0 ){
        inside <- fields::in.poly.grid(gridList, polyStd)
        if (sum(inside) > 0) {
          # note basisScale factor has to be added because this sum is
          # over the  standard basis functions with  scale  1.0
          tmpIntegral <-sum(zGrid[inside])  * dx * dy * basisScale ^ 2
          # accumulate sparse matrix information
          #J <- c(J, j)
          #K <- c(K, k + basisOffset) # offset adjusts for preceding levels. 
          tempIntegral <- c(tempIntegral, tmpIntegral) 
        }
      }
      
    }
    tempIntegral
  }
}

parallel::stopCluster(cl = my.cluster)


#test section

#base for loop test how to construct the correct data object
for (j in 1:N1) {
  # loop over the polygon regions (This could be parallelized using a foreach constuct. )
  polyTmp <- (polyGroups)[[j]]
  
  for (L in 1:nLevel) {
    # get info for the L^th level of the multi-resolution basis. 
    basisIndex <-
      (1:LKinfo$latticeInfo$mLevel[L])
    basisOffset<- LKinfo$latticeInfo$offset[L]
    basisScale <- LKinfo$latticeInfo$delta[L]*
      LKinfo$basisInfo$overlap
    basisGridList <- LKinfo$latticeInfo$grid[[L]]
    basisCenters <- make.surface.grid(basisGridList)
    
    for (k in basisIndex){
      # k is position at the L^th level
      polyStd <-
        cbind((polyTmp[, 1] - basisCenters[k, 1]) / basisScale,
              (polyTmp[, 2] - basisCenters[k, 2]) / basisScale)
      # points in the region
      
      allOutside<- in.poly(  polyStd, boundingBox)
      # only look at polygonal if at least one point is inside the bounding box
      # note that integrals with basis function support outside of the
      # polygon are skipped and by defualt set to zero
      
      if(sum(allOutside) > 0 ){
        inside <- in.poly.grid(gridList, polyStd)
        if (sum(inside) > 0) {
          # note basisScale factor has to be added because this sum is
          # over the  standard basis functions with  scale  1.0
          tmpIntegral <- sum(zGrid[inside])  * dx * dy * basisScale ^ 2
          # accumulate sparse matrix information
          J <- c(J, j)
          K <- c(K, k + basisOffset) # offset adjusts for preceding levels. 
          integral <- c(integral, tmpIntegral)
        }
      }
      
    }
  }
}
    
obj<- list( ind= cbind( J,K), ra= integral, 
                da= c( N1,N2))

X <- spind2spam(obj)

#---test-----###

#create empty df with columns for integrals, J, and K (where J and K are used in sparse)
basisInt <- as.data.frame(matrix(ncol = 3))

#cluster set-up
n.cores <- cores

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)

#TODO: fix the the combine method and correct output
basisInt <- foreach(j = 1:N1, .combine=rbind) %dopar% {
  polyTmp <- polyGroups[[j]]
  
  integrals <- as.data.frame(matrix(ncol = 3))
  for (l in 1:nLevel) {
    # get info for the L^th level of the multi-resolution basis. 
    basisIndex <-(1:LKinfo$latticeInfo$mLevel[l])
    basisOffset <- LKinfo$latticeInfo$offset[l]
    basisScale <- LKinfo$latticeInfo$delta[l]*LKinfo$basisInfo$overlap
    basisGridList <- LKinfo$latticeInfo$grid[[l]]
    basisCenters <- fields::make.surface.grid(basisGridList) #requires fields
    
    for (k in basisIndex){
      # k is position at the L^th level
      polyStd <-
        cbind((polyTmp[ ,1] - basisCenters[k, 1]) / basisScale,
              (polyTmp[ ,2] - basisCenters[k, 2]) / basisScale)
      # points in the region
      
      allOutside<- fields::in.poly(  polyStd, boundingBox)
      # only look at polygonal if at least one point is inside the bounding box
      # note that integrals with basis function support outside of the
      # polygon are skipped and by defualt set to zero
      
      if(sum(allOutside) > 0 ){
        inside <- fields::in.poly.grid(gridList, polyStd)
        if (sum(inside) > 0) {
          # note basisScale factor has to be added because this sum is
          # over the  standard basis functions with  scale  1.0
          tmpIntegral <- sum(zGrid[inside])  * dx * dy * basisScale ^ 2
          # accumulate sparse matrix information
          tempbasis_df <- c(tmpIntegral, j,  k + basisOffset)
          # offset adjusts for preceding levels. 
          integrals <- rbind(integrals, tempbasis_df)
          }
        }
    
    }
  }
  integrals[-1, ]
}

parallel::stopCluster(cl = my.cluster)

colnames(basisInt) <- c("int", "J", "K")

#now are they the same as for loop?
int2 <- basisInt$int
J2 <- basisInt$J
K2 <- basisInt$K

test.for.zero(int2, integral)
test.for.zero(J2, J)
test.for.zero(K2, K)


###-------------####
#test with new function:


#set up Lkinfo: (normalize must be false)
sDomain <- sDomain<- rbind( c( -1,-1),
                            c( 1,2))

#look up alpha (setting nu as 1 instead)
LKinfo <- LKrigSetup(sDomain, NC=25, nlevel=2, a.wght = 4.1,
                     NC.buffer = 0, normalize = FALSE, nu = 1)
#check alpha setup for multiple levels, missing alpha specification
##see paper for more details on alpha selection.

M = 200

source("DF_LK/R/basisIntegralNew.R")
source("DF_LK/R/basisIntegral.R")

library(tictoc)

tic()
X <- basisIntegral(polyGroups, LKinfo, M = 200)
toc()

tic()
X <- basisIntegralForEach(polyGroups, LKinfo, M = 200, cores= 8)
toc()
