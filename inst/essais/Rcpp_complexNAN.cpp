#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::ComplexVector test() {
  Rcomplex z;
  z.r = NAN;
  z.i = 0;
  return Rcpp::ComplexVector::create(z);
}

// [[Rcpp::export]]
Rcpp::ComplexVector zero() {
  Rcomplex z; // ok but warning
  return Rcpp::ComplexVector::create(z);
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
test()
*/
