#include "RcppArmadillo.h"

using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


//' Perform the forward backward algorithm (see Rabiner 89)
//'
//' @param m the number of positions (hypothesis)
//' @param A a matrix 2 * 2 the transition probabilities
//' @param f0x a vector of the values of the density under the null hypothesis on the observations
//' @param f1x a vector of the values of the density under the alternative hypothesis on the observations
//' @param Pi a vector of the initial state probabilities
//'
//' @return alpha the forward variables, the lines corespond to the position, the first column is for state 0 and the second for state one
//' @return beta the backward variables (same as alpha)
//' @return gamma  matrix such that gamma\[i , 0\] = \eqn{P(\theta_i = 0 | X)}, gamma\[i,1\]  = \eqn{P(\theta_i = 1 | X)}
//' @return ksi  matrix such that ksi\[i , 0\] = \eqn{P(\theta_i = 0 | X)}, ksi\[i,1\]  = \eqn{P(\theta_i = 1 | X)}
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
// [[Rcpp::export]]
arma::vec  viterbi(int m,
                    arma::mat A,
                    arma::vec f0x,
                    arma::vec f1x,
                    arma::vec Pi) {
  arma::mat Delta(m, 2);
  arma::mat Psi(m, 2);

  Delta(0, 0) = Pi(0) * f0x(0);
  Delta(0, 1) = Pi(1) * f1x(0);
  Psi(0, 0)   = 0;
  Psi(0, 1)   = 0;
  for (int i = 1; i < m ; ++i)
  {
    vec value0(2);
    vec value1(2);
    value0(0) = Delta(i - 1, 0) * A(0, 0) * f0x(i);
    value0(1) = Delta(i - 1, 1) * A(1, 0) * f0x(i);
    value1(0) = Delta(i - 1, 0) * A(0, 1) * f1x(i);
    value1(1) = Delta(i - 1, 1) * A(1, 1) * f1x(i);
    Psi(i, 0) = value0.index_max();
    Psi(i, 1) = value1.index_max();
    Delta(i, 0) = max(value0);
    Delta(i, 1) = max(value1);
  }

  arma::vec best_path(m);
  best_path(m - 1) = Delta.row(m - 1).index_max();
//
  for(int i = m - 2; i > - 1; --i){
    best_path(i) = Psi(i + 1, best_path(i + 1));
  }


  return best_path;
}


//' Perform the forward backward algorithm (see Rabiner 89)
//'
//' @param m the number of positions (hypothesis)
//' @param A a matrix 2 * 2 the transition probabilities
//' @param f0x a vector of the values of the density under the null hypothesis on the observations
//' @param f1x a vector of the values of the density under the alternative hypothesis on the observations
//' @param Pi a vector of the initial state probabilities
//'
//' @return alpha the forward variables, the lines corespond to the position, the first column is for state 0 and the second for state one
//' @return beta the backward variables (same as alpha)
//' @return gamma  matrix such that gamma\[i , 0\] = \eqn{P(\theta_i = 0 | X)}, gamma\[i,1\]  = \eqn{P(\theta_i = 1 | X)}
//' @return ksi  matrix such that ksi\[i , 0\] = \eqn{P(\theta_i = 0 | X)}, ksi\[i,1\]  = \eqn{P(\theta_i = 1 | X)}
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
// [[Rcpp::export]]
arma::vec  viterbi_log(int m,
                   arma::mat A_log,
                   arma::vec f0x_log,
                   arma::vec f1x_log,
                   arma::vec Pi_log) {
  arma::mat Phi(m, 2);
  arma::mat Psi(m, 2);

  Phi(0, 0) =   Pi_log(0) + f0x_log(0);
  Phi(0, 1) = Pi_log(1) + f1x_log(0);
  Psi(0, 0)   = 0;
  Psi(0, 1)   = 0;
  for (int i = 1; i < m ; ++i)
  {
    vec value0(2);
    vec value1(2);
    value0(0) =  Phi(i - 1, 0) + A_log(0, 0) + f0x_log(i);
    value0(1) =  Phi(i - 1, 1) + A_log(1, 0) + f0x_log(i);
    value1(0) =  Phi(i - 1, 0) + A_log(0, 1) + f1x_log(i);
    value1(1) =  Phi(i - 1, 1) + A_log(1, 1) + f1x_log(i);
    Psi(i, 0) = value0.index_max();
    Psi(i, 1) = value1.index_max();
    Phi(i, 0) = max(value0);
    Phi(i, 1) = max(value1);
  }

  arma::vec best_path(m);
  best_path(m - 1) = Phi.row(m - 1).index_max();
  for(int i = m - 2; i > - 1; --i){
    best_path(i) = Psi(i + 1, best_path(i + 1));
  }


  return best_path;
}

