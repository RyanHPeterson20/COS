
#bspline polynomial function

#TODO: add in some checks for outside of support so we don't have to evaluate everything
#TODO: convert to C (check which runs faster)

#compare pmax and max
#cubicSplineBasis <- function(x){
#  return((1/6) *(pmax(x + 2, 0)^3 - 4 * (pmax(x + 1, 0)^3) + 6 * (pmax(x, 0)^3) -
#                   4 * (pmax(x - 1, 0)^3) + pmax(x - 2, 0)^3))  
#}

#testing variation and new name
CubicBSpline <- function(d){
  d <- d*2 #added in to deal with 
  return((1/6) *(pmax(d + 2, 0)^3 - 4 * (pmax(d + 1, 0)^3) + 6 * (pmax(d, 0)^3) -
                   4 * (pmax(d - 1, 0)^3) + pmax(d - 2, 0)^3) * (abs(d) < 2))  
}


