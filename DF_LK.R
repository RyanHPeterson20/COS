#script for running data fusion with LatticeKrig

##TODO list:
#-- Create 

suppressMessages(library(LatticeKrig))

#parallel librarues
suppressMessages(library(foreach))
suppressMessages(library(parallel))
suppressMessages(library(doParallel))


#load scripts for functions
source("COSExample/R/basisIntegral.R")
source("COSExample/R/integralPolyFunction.R")
# patch in newer version of LK functions
source("COSExample/R/LatticeKrig.R")
source("COSExample/R/LKrigFindLambda.R")
source("COSExample/R/print.LatticeKrig.R" )
# patch in updates to integral functions
source("DF_LK/R/basisIntegralNew.R")
source("DF_LK/R/integralPolyFunctionNew.R")

#load data 
load("DF_LK/DFSynth_RandomBig.rda")



#updated domain
sDomain <- sDomain<- rbind( c( -1,-1),
                            c( 1,2))

#set up Lkinfo: (normalize must be false)
#look up alpha (setting nu as 1 instead)
LKinfo <- LKrigSetup(sDomain, NC=25, nlevel=2, a.wght = 4.1,
                     NC.buffer = 2, normalize = FALSE, nu = 1)

#creates a couple functions
FUNX<- function(s){
  s[,1]
}

FUNY<- function(s){
  s[,2]
}

#------------------------------------------#

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
X_1<- basisIntegralForEach( newpolyGroupsBig, LKinfo, M = 400, cores = 8)
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

## for data set 2

U1<- integralPolyFunctionForEach(newpolyGroupsBig_2, M = 200, cores = 8)  # default is constant function 
U2<- integralPolyFunctionForEach(newpolyGroupsBig_2,
                          FUN=function(s){s[,1]}, M = 200, cores =8)
U3<- integralPolyFunctionForEach(newpolyGroupsBig_2,
                          FUN=function(s){s[,2]}, M = 200, cores = 8)
U_2<- cbind( U1,U2,U3) # U can be a dense matrix


# integrals of basis functions
X_2<- basisIntegralForEach( newpolyGroupsBig_2, LKinfo, M = 200, cores = 8)
X_2<- spind2spam(X_2)

# estimating lambda and the a.wght parameter
fit2<- LatticeKrig( sDomain, newobsRandomBig_2, U=U_2, X=X_2, LKinfo=LKinfo,
                    findAwght=TRUE) #does Lkinfo overwrite findAweight?

set.panel(1,2)
zlim<- range( gTrue)
surface(as.surface( gridTrue, gTrue) , 
        zlim =zlim, col=viridis(256))
title("true surface")
surface( fit2, zlim =zlim)
title("lambda and a.wght MLEs")

simOut2<- LKrig.sim.conditional( fit2,  M=100) 

set.panel(1,2)
zlim<- range( gTrue)
surface( fit2, col=turbo(256))
title("lambda and a.wght MLEs")
imagePlot(as.surface(simOut2$x.grid,simOut2$SE))
title("Prediction SE")



##FUSION###

#Base DF with LatticeKrig
U_F <- rbind(U_1, U_2)

#create new matrix X_F
X_F <- rbind(X_1, X_2)


y_full <- append(newobsRandomBig, newobsRandomBig_2)

# estimating lambda and the a.wght parameter
fit_F<- LatticeKrig( sDomain, y_full, U=U_F, X=X_F, LKinfo=LKinfo,
                     findAwght=TRUE)

set.panel(1,2)
zlim<- range( gTrue)
surface(as.surface( gridTrue, gTrue), col=viridis(256))
title("true surface")
surface( fit_F)
title("lambda and a.wght MLEs")

simOut_F<- LKrig.sim.conditional( fit_F,  M=100) 

set.panel(1,2)
zlim<- range( gTrue)
surface( fit_F, zlim = zlim, col=turbo(256))
title("lambda and a.wght MLEs")
imagePlot(as.surface(simOut_F$x.grid,simOut_F$SE))
title("Prediction SE")

save(U_1, X_1,
     U_2, X_1,
     U_F, X_F,
     file = "DF_matrices.rda")

save(fit1, simOut1,
     fit2, simOut2,
     fit_F, simOut_F,
     file = "DF_predictions.rda")