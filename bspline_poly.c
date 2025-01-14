#include <R.h>
#include <Rinternals.h>
#include <math.h>

// comment 
// testing power functions
void cubicSplineBasis(double *x, double *result){
  
  double term1 = pow(fmax(*x + 2, 0), 3);
  double term2 = 4 * pow(fmax(*x + 1, 0), 3);
  double term3 = 6 * pow(fmax(*x, 0), 3);
  double term4 = 4 * pow(fmax(*x - 1, 0), 3);
  double term5 = pow(fmax(*x - 2, 0), 3);
  
  *result = (1.0 / 6.0) * (term1 - term2 + term3 - term4 + term5);
}

