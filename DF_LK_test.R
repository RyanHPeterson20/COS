#test section for new methods and additions to the COS LK code
#uses simplified block data (later we'll include more complex data to test.)

#libraries
suppressMessages(library(LatticeKrig))

#parallel librarues
suppressMessages(library(foreach))
suppressMessages(library(parallel))
suppressMessages(library(doParallel))


#load scripts for functions
setwd("~/COS_LK")

source("COSExample/R/basisIntegral.R")
source("COSExample/R/integralPolyFunction.R")
# patch in newer version of LK functions
source("COSExample/R/LatticeKrig.R")
source("COSExample/R/LKrigFindLambda.R")
source("COSExample/R/print.LatticeKrig.R" )
# patch in updates to integral functions
source("DF_LK/R/basisIntegralNew.R")
source("DF_LK/R/integralPolyFunctionNew.R")

#load test data:
load("DF_LK/DFSynth_RandomBig.rda")
load("DF_LK/DFSynthData_Base.rda")

#updated domain
sDomain <- sDomain<- rbind( c( -1,-1),
                            c( 1, 2))

#set up Lkinfo: (normalize must be false)
#look up alpha (setting nu as 1 instead)
LKinfo <- LKrigSetup(sDomain, NC=15, nlevel=1, a.wght = 4.1,
                     NC.buffer = 0, normalize = FALSE, nu = 1)

#creates a couple functions
FUNX<- function(s){
  s[,1]
}

FUNY<- function(s){
  s[,2]
}


#perform the COSP LK for each data set, and get SE surface as well

#begin with regular tiles for predictions, we'll let LK find lambda and aRange
# create matrix for fixed effects integrals
# (a linear function of coordinates in this case)
U1<- integralPolyFunctionForEach(newpolyGroupsBig, M = 200, cores = 8)  # default is constant function 
U2<- integralPolyFunctionForEach(newpolyGroupsBig,
                                 FUN=function(s){s[,1]}, M = 200, cores = 8)
U3<- integralPolyFunctionForEach(newpolyGroupsBig,
                                 FUN=function(s){s[,2]}, M = 200, cores = 8)
U_1<- cbind( U1,U2,U3) # U can be a dense matrix


# integrals of basis functions
X_1<- basisIntegral_New( newpolyGroupsBig, LKinfo, M = 200, cores = 8, normalize = TRUE)
X_1<- spind2spam(X_1)

# estimating lambda and the a.wght parameter
fit1<- LatticeKrig( sDomain, newobsRandomBig, U=U_1, X=X_1, LKinfo=LKinfo,
                    findAwght=TRUE)

set.panel(1,2)
zlim<- range( gTrue)
surface(as.surface( gridTrue, gTrue) , 
        zlim =zlim, col=viridis(256))
title("true surface")
surface( fit1, zlim =zlim)
title("lambda and a.wght MLEs")


simOut1<- LKrig.sim.conditional( fit1,  M=100) 

set.panel(1,2)
zlim<- range( gTrue)
surface( fit1, zlim=zlim)
title("lambda and a.wght MLEs")
imagePlot(as.surface(simOut1$x.grid,simOut1$SE))
title("Prediction SE")


# integrals of basis functions
X_1<- basisIntegral_New( newpolyGroupsBig, LKinfo, M = 200, cores = 8, normalize = FALSE)
X_1<- spind2spam(X_1)

# estimating lambda and the a.wght parameter
fit1<- LatticeKrig( sDomain, newobsRandomBig, U=U_1, X=X_1, LKinfo=LKinfo,
                    findAwght=TRUE)

set.panel(1,2)
zlim<- range( gTrue)
surface(as.surface( gridTrue, gTrue) , 
        zlim =zlim, col=viridis(256))
title("true surface")
surface( fit1, zlim =zlim)
title("lambda and a.wght MLEs")


simOut1<- LKrig.sim.conditional( fit1,  M=100) 

set.panel(1,2)
zlim<- range( gTrue)
surface( fit1, zlim=zlim)
title("Normalize = FALSE")
imagePlot(as.surface(simOut1$x.grid,simOut1$SE))
title("Prediction SE")



