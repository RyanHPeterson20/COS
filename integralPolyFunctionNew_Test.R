#delete library text/test later
#both required
library(foreach)
library(parallel)

detectCores()

# parameters:
## polyGroups
## FUN (default NULL) integration function
## M (default 200) Integration points
## 

#requires: foreach

#returns: 

forloop_action <- function(polyGroups, j, FUN, M){
  
  
  polyTmp <- (polyGroups)[[j]]
  polyGridList<- list(
    x= seq(min(polyTmp[,1]), max(polyTmp[,1]), length.out=M),
    y= seq(min(polyTmp[,2]), max(polyTmp[,2]), length.out=M)
  )
  dxP<-  polyGridList$x[2] - polyGridList$x[1]
  dyP<-  polyGridList$y[2] - polyGridList$y[1]
  
  polyGrid<- make.surface.grid( polyGridList)
  ind<- in.poly(polyGrid, polyTmp)
  polyGrid<- polyGrid[ind,]
  if( !is.null(FUN)){
    look <- FUN(polyGrid, ... )
    theIntegrals[j]<- sum( look)*dxP*dyP
  }
  else{
    # just find the area
    theIntegrals[j]<- sum( ind)*dxP*dyP
  }
  return()
}

integralPolyFunction <- function(polyGroups,
                                 FUN=NULL,
                                 M = 200,...) {
  
  N1 <- length(polyGroups)
  theIntegrals<-rep(NA, N1)
  #TODO: add in for each
  foreach (j = 1:N1) %do% forloop_action(polyGroups, j, FUN)
    
  return(theIntegrals)
}


#test block----------------------
foreach(i=1:3) %do% sqrt(i)

m <- matrix(rnorm(9), 3, 3)
foreach(i=1:ncol(m), .combine=c) %do%
  mean(m[,i])

#try parallel
foreach(i=1:nrow(m), .combine=rbind) %dopar%
  (m[i,] / mean(m[i,]))

#---------------------------------
load("DF_LK/COSExampleNew.rda")

#use polyGroups2 for testing
N1 <- length(polyGroups2)
M <- 100
FUN <- 

test_function <- function(polyGroups2, M){
  polyTmp <- polyGroups2
  
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
  
  # just find the area
  integral <- sum( ind)*dxP*dyP

  return(integral)
}

N1 <- length(polyGroupsRandomBig)

num_cores <- detectCores()

# Create a cluster with multiple workers
cl <- makeCluster(num_cores)

# Check the number of workers in the cluster
num_workers <- length(clusterApply(cl, 1, function(x) Sys.info()[['nodename']]))
cat("Number of workers in use:", num_workers, "\n")



system.time(
theIntegrals <- foreach(i=1:N1, M = rep(200, N1), .combine=rbind) %dopar% test_function(polyGroupsRandomBig[[i]], M)
)

# Close the cluster
stopCluster(cl)


forloop_action(polyGroups2[[5]], M = 200, FUN=function(s){s[,1]})


FUN=function(s){s[,1]}
FUN=function(s){s[,2]}
look <- FUN(polyGrid)
sum( look)*dxP*dyP


#use this viz later on,
## to visualize the integration points for a single block.
## use a similar viz for overlap
plot(polyGrid, pch = 16, cex = 0.5)



#test new integralPolyFunction----------------

library(LatticeKrig)
library(fields)
library(foreach)
library(parallel)


source("DF_LK/R/integralPolyFunctionNew.R")
source("DF_LK/R/integralPolyFunction.R")

load("DF_LK/COSExampleNew.rda")

#create covariate matrix from 
FUNX<- function(s){
  s[,1]
}

FUNY<- function(s){
  s[,2]
}

system.time(
  U1<- integralPolyFunction(polyGroupsRandomBig, M = 200)  # default is constant function 
)

system.time(
  U1<- integralPolyFunctionNew(polyGroupsRandomBig, M = 200)  # default is constant function 
)

U2<- integralPolyFunction(polyGroupsRandomBig,
                          FUN=function(s){s[,1]}, M = 200)
U3<- integralPolyFunction(polyGroupsRandomBig,
                          FUN=function(s){s[,2]}, M = 200)
U_1<- cbind( U1,U2,U3)

system.time(
U1<- integralPolyFunctionNew(polyGroupsRandomBig, M = 200)  # default is constant function 
)

U2<- integralPolyFunctionNew(polyGroupsRandomBig,
                          FUN=function(s){s[,1]}, M = 200)
U3<- integralPolyFunctionNew(polyGroupsRandomBig,
                          FUN=function(s){s[,2]}, M = 200)
U_1<- cbind( U1,U2,U3)



