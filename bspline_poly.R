
#bspline polynomial function

cubicSplineBasis <- function(x){
  return((1/6) *(pmax(x + 2, 0)^3 - 4 * (pmax(x + 1, 0)^3) + 6 * (pmax(x, 0)^3) -
                   4 * (pmax(x - 1, 0)^3) + pmax(x - 2, 0)^3))  
}
