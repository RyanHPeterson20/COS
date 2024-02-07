
integralPolyFunctionNew <- function(polyGroups,
                                 FUN=NULL,
                                 M = 200,...) {
  # setup basis function info
  
  #TODO: add in library check error catch
  
  forEach_function <- function(polyGroups, M, FUN){
    polyTmp <- polyGroups
    
    polyGridList<- list(
      x= seq(min(polyTmp[,1]), max(polyTmp[,1]), length.out=M),
      y= seq(min(polyTmp[,2]), max(polyTmp[,2]), length.out=M)
    )
    
    #integral delta
    dxP<-  polyGridList$x[2] - polyGridList$x[1]
    dyP<-  polyGridList$y[2] - polyGridList$y[1]
    
    polyGrid<- make.surface.grid( polyGridList)
    ind <- in.poly(polyGrid, polyTmp)
    polyGrid<- polyGrid[ind,]
    
    
    if( !is.null(FUN)){
      look <- FUN(polyGrid, ... )
      inntegral <- sum( look)*dxP*dyP
    }
    else{
      # just find the area
      integral <- sum( ind)*dxP*dyP
    }
    return(integral)
  }

  N1 <- length(polyGroups)
  
  theIntegrals <- foreach(i=1:N1, M = rep(M, N1), .combine=rbind) %dopar% 
    forEach_function(polyGroups[[i]], M, FUN=FUN)
  
  return(theIntegrals)
}  