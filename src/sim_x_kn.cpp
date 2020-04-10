#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' Simulate a processus from  an heterogeneuous markov chain. This heterogeneous markov Chain is the low of the state of an HMM given the obervation.
//'
//' @param alpha a matrix m * 2  containing the forward variables
//' @param beta a matrix m * 2  containing the backward variables
//' @param A a matrix 2 * 2 the transition probabilities
//' @param f0x a vector of the values of the density under the null hypothesis on the observations
//' @param f1x a vector of the values of the density under the alternative hypothesis on the observations
//' @param Pi vector of the initial state probabilities for each state.
//'
//' @return Product of matrices
//' @export
//'
//' @examples
//'  m <-  10
//'  A <- matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T)
//'  f0 <- c(0, 1)
//'  f1 <- c(2, 1)
//'  Pi <- c( 0.9, 0.1)
//'  rdata <- simulate.data.hmm.2states(m, Pi, A, f0, f1)
//'  x <- rdata$x
//'  theta <- rdata$theta
//'  mod <- for_back(m, A, f0x, f1x, Pi)
//'  f0x <- dnorm(x, f0[1], f0[2])
//'  f1x <- dnorm(x, f1[1], f1[2])
//'  alpha <- mod$alpha
//'  beta <- mod$beta
//'  sim_x_kn(m, alpha, beta, A, Pi, f0x, f1x)
// [[Rcpp::export]]
NumericVector sim_x_kn(int m, NumericMatrix alpha,
                       NumericMatrix beta,
                       NumericMatrix A,
                       NumericVector Pi,
                       NumericVector f0x,
                       NumericVector f1x) {
  NumericVector out(m);
  double ai12;
  double ai22;
  double li0;
  double li1;
  out[0] = as<double>(rbinom(1, 1, Pi[1]));
  double un = 1;
  for(int i = 1; i < m; ++i) {
    li0 = alpha(i - 1 , 0) * beta(i - 1 , 0) / (alpha(i - 1 , 0) * beta(i - 1 , 0) + alpha(i - 1 , 1) * beta(i - 1 , 1));
    li1 = alpha(i - 1 , 1) * beta(i - 1 , 1) / (alpha(i - 1 , 0) * beta(i - 1 , 0) + alpha(i - 1 , 1) * beta(i - 1 , 1));
    ai12 =  beta(i, 1) * alpha(i - 1, 0) * f1x[i] * A(0, 1) / li0;
    ai22 =  beta(i, 1) * alpha(i - 1, 1) * f1x[i] * A(1, 1) / li1;
    out[i] =  (1 - out[i - 1]) * as<double>(rbinom(1, 1, std::min<double>(ai12,1.0))) +(out[i - 1]) * as<double>(rbinom(1, 1, std::min<double>(ai22,1.0)));
  }
  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
# timesTwo(42)
*/
