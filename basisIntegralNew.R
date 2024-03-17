# parameters:
## polyGroups
## LKinfo
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
  zGrid <- Wendland(sqrt(xyGrid[, 1] ^ 2 + xyGrid[, 2] ^ 2),
                    dimension = 2,
                    k = 2)
  stdIntegral<- sum( zGrid)*dx*dy
  
  
  ##section pending completion in test section:
  
  
  ##integration section: in parallel
  
  #cluster set-up
  n.cores <- cores
  
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  doParallel::registerDoParallel(cl = my.cluster)
  
  basisInt <- foreach(j = 1:N1, .combine=rbind) %dopar% {
    polyTmp <- polyGroups[[j]]
  }
  
  parallel::stopCluster(cl = my.cluster)
  
  return(basisInt)
  
}  