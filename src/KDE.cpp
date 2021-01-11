#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector lapin(NumericVector x, Function f) {
  NumericVector res = f(x);
  return res;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
lapin(10, dnorm)
*/
