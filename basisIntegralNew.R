#required libraries: 
#suppressMessages(library(foreach))
#suppressMessages(library(parallel))
#suppressMessages(library(doParallel))

# parameters:
## polyGroups
## LKinfo
## Cores
## M (default 200) Integration points
## check (default TRUE)

basisIntegralForEach <- function(polyGroups,
                          LKinfo,
                          M = 200, cores, check=TRUE) {
  
  # setup basis function info
  
  nLevel <- LKinfo$nlevel
  
  gridList <- list(x = seq(-1, 1, length.out = M),
                   y = seq(-1, 1, length.out = M))
  
  #TODO: give explanation of this
  boundingBox<- rbind( c( -1,-1),
                       c( 1,-1),
                       c( 1,1),
                       c( -1, 1),
                       c( -1,-1))
  
  dx<- gridList$x[2]- gridList$x[1]
  dy<- gridList$y[2]- gridList$y[1]
  
  
  # integral for standard basis function. (is this necessary?)
  xyGrid <- make.surface.grid(gridList)
  
  # values for the standard basis function on the grid
  
  #TODO: replace this with an option for radial or tensor
  ## evaluate  Radial.basis() or Tensor.basis() at some grid centered on (0,0)
  #pull basystype fron lkinfo
  #make sure this is pretty standardized
  ## see LKrig.basis for code example
  if (condition) {
   #radial 
  }
  if (condition) {
   #tensor
  }
  
  zGrid <- Wendland(sqrt(xyGrid[, 1] ^ 2 + xyGrid[, 2] ^ 2),
                    dimension = 2,
                    k = 2)
  stdIntegral<- sum( zGrid)*dx*dy
  
  
  ##section pending completion in test section:
  #
  N1 <- length(polyGroups)
  # total number of basis functions
  N2 <- LKinfo$latticeInfo$m
  # basis function scales
  delta <- LKinfo$latticeInfo$delta
  
  checkIntegral<- NULL # stays NULL if check == FALSE
  
  ##integration section: in parallel
  
  #create empty df with columns for integrals, J, and K (where J and K are used in sparse)
  basisInt <- as.data.frame(matrix(ncol = 3))
  
  #cluster set-up
  n.cores <- cores
  
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  doParallel::registerDoParallel(cl = my.cluster)
  
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
  
  # spind sparse format.
  obj<- list( ind= cbind( basisInt$J,basisInt$K), ra= basisInt$int, 
              da= c( N1,N2))
  
  return(obj)
  
}  