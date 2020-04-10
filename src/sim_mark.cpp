#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' Simple matrix multiplication
//'
//' @param A Matrix
//' @param B Matrix
//'
//' @return Product of matrices
//' @export
//'
//' @examples
//' A <- matrix(1:9, 3, 3)
//' B <- matrix(11:19, 3, 3)
//' matrix_mult_cpp(A, B)
// [[Rcpp::export]]
arma::vec sim_markov(int m, arma::vec Pi, arma::mat A) {
 arma::vec theta(m);
  theta(0) = as<double>(rbinom(1, 1, Pi(1)));
  for(int i = 1; i < m; ++i) {
    theta(i) =  (1 - theta(i - 1)) * as<double>(rbinom(1, 1, A(0, 1))) +
      theta(i - 1) * as<double>(rbinom(1, 1, A(1, 1)));
  }
  }
