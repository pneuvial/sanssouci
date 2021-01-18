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
//'  rdata <- sim_hmm_2states(m, Pi, A, f0, f1)
//'  x <- rdata$x
//'  theta <- rdata$theta
//'  mod <- for_back(m, A, f0x = dnorm(x), f1x = dnorm(x, 1,2), Pi)
// [[Rcpp::export]]
Rcpp::List for_back(int m,
                    arma::mat A,
                    arma::vec f0x,
                    arma::vec f1x,
                    arma::vec Pi) {
  arma::mat alpha(m, 2);
  arma::vec c0(m);

  c0(0) = 1 / (Pi(0)*f0x(0) + Pi(1)*f1x(0));
  alpha(0, 0) = Pi(0)*f0x(0) * c0(0);
  alpha(0, 1) = Pi(1)*f1x(0) * c0(0);


  for (int i = 1; i < m ; ++i)
  {

    alpha(i, 0)  = f0x(i) * (alpha(i - 1, 0) * A(0, 0) + alpha(i - 1, 1) * A(1, 0));
    alpha(i, 1)  = f1x(i) * (alpha(i - 1, 0) * A(0, 1) + alpha(i -1, 1) * A(1, 1));
    c0(i) = 1 / ( alpha(i, 0) +  alpha(i, 1));
    alpha(i, 0) = alpha(i, 0) * c0(i);
    alpha(i, 1) = alpha(i, 1) * c0(i);
  }
  arma::mat beta(m, 2);
  beta(m - 1, 0) = 1;
  beta(m - 1, 1) = 1;

  beta.row(m - 1) = c0(m - 1) * beta.row(m - 1);
  arma::mat gamma(m, 2);
  gamma(m - 1, 0) = alpha(m - 1, 0) * beta(m - 1 ,0) / ( alpha(m - 1, 0) * beta(m - 1 ,0) + alpha(m - 1, 1) * beta(m - 1 ,1));
  gamma(m - 1, 1) = alpha(m - 1, 1) * beta(m - 1 ,1) / ( alpha(m - 1, 0) * beta(m - 1 ,0) + alpha(m - 1, 1) * beta(m - 1 ,1));
  arma::mat ksi(m, 4);
  for (int i = m - 2; i > -1; --i)
  {
    beta (i, 0) = c0(i) * (A(0, 0) * f0x(i + 1) * beta(i + 1, 0) + A(0, 1) * f1x(i + 1) * beta(i + 1, 1));
    beta (i, 1) = c0(i) * (A(1, 0) * f0x(i + 1) * beta(i + 1, 0) + A(1, 1) * f1x(i + 1) * beta(i + 1, 1));
    gamma(i, 0) = alpha(i, 0) * beta(i ,0) / ( alpha(i, 0) * beta(i ,0) + alpha(i, 1) * beta(i ,1));
    gamma(i, 1) = alpha(i, 1) * beta(i ,1) / ( alpha(i, 0) * beta(i ,0) + alpha(i, 1) * beta(i ,1));
    double denom =0;
    for(int j = 0; j < 2; j++){
      for(int k = 0; k < 2; k++){
        double fx = (1-k)*f0x(i+1)+k*f1x(i+1);
        denom =    alpha(i, j) * A(j, k) * fx * beta(i + 1, k) + denom;
      }
    }

    for(int j = 0; j < 2; j++){
      ksi(i, j) = alpha(i, j) * A(j, 0) * f0x(i + 1) * beta(i + 1, 0) /  denom;
      ksi(i, j  + 2) = alpha(i, j) * A(j, 1) * f1x(i + 1) * beta(i + 1, 1) / denom;
    }
  }


  return Rcpp::List::create(Rcpp::Named("alpha")=alpha,
                            Rcpp::Named("beta")= beta,
                            Rcpp::Named("gamma")= gamma,
                            Rcpp::Named("ksi")= ksi);
}




//' Use EM algorithm to estimate the parameters A and Pi of an HMM
//'
//' @param m the number of positions (hypothesis)
//' @param alpha a matrix m * 2  containing the forward variables
//' @param beta a matrix m * 2  containing the backward variables
//' @param A a matrix 2 * 2 the transition probabilities
//' @param f0x a vector of the values of the density under the null hypothesis on the observations
//' @param f1x a vector of the values of the density under the alternative hypothesis on the observations
//' @param Pi a vector of the initial state probabilities
//' @param eps the value ta reach for the convergence
//' @param maxit integer, the maximum number of iteration
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
// [[Rcpp::export]]
Rcpp::List Em_hmm(int m,
                  arma::mat A,
                  arma::vec f0x,
                  arma::vec f1x,
                  arma::vec Pi,
                  double eps,
                  int maxit) {
  arma::vec diff = vec(2);
  diff(0) = eps + 1;
  int i = 0;
  while((max(diff) > eps & i < maxit) )
  {
    arma::mat A_old = A;
    arma::vec Pi_old = Pi;
    Rcpp::List b_f = for_back(m, A, f0x, f1x, Pi);
    arma::mat ksi = b_f["ksi"];
    arma::mat gamma = b_f["gamma"];
    arma::rowvec sum_gamma = sum(gamma.rows(0,m-2));
    arma::rowvec sum_ksi = sum(ksi.rows(0,m-2));
    Pi = trans(gamma.row(1));
    for(int j = 0; j < 2; j++){
      for(int k = 0; k < 2; k++){
        A(j, k) = sum_ksi(j + 2 * k) / sum_gamma(j);
      }
    }
    diff(0) = abs(A - A_old).max();
    diff(1) = max(Pi - Pi_old);
    i++;
  }
  return Rcpp::List::create(Rcpp::Named("A") = A,
                            Rcpp::Named("Pi") = Pi,
                            Rcpp::Named("i") = i);
}




//' Use EM algorithm to estimate the parameters A and Pi of an HMM
//'
//' @param m the number of positions (hypothesis)
//' @param alpha a matrix m * 2  containing the forward variables
//' @param beta a matrix m * 2  containing the backward variables
//' @param A a matrix 2 * 2 the transition probabilities
//' @param f0x a vector of the values of the density under the null hypothesis on the observations
//' @param f1x a vector of the values of the density under the alternative hypothesis on the observations
//' @param Pi a vector of the initial state probabilities
//' @param eps the value ta reach for the convergence
//' @param maxit integer, the maximum number of iteration
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
// [[Rcpp::export]]
Rcpp::List Em_f1(int m,
                  arma::mat A,
                  arma::vec Pi,
                  arma::vec f0x,
                  arma::vec f1x,
                  Rcpp::List fw_bc_EM,
                  arma::vec x,
                  double eps,
                  int maxit,
                  double h) {
  arma::vec diff = vec(3);
  diff(0) = eps + 1;
  int i = 0;
  Rcpp::List b_f;
  while(((max(diff) > eps) & (i < maxit)))
  {
    arma::mat A_old = A;
    arma::vec Pi_old = Pi;
    arma::vec f1x_old = f1x;
    Rcpp::List  EM = Em_hmm(m , A, f0x, f1x,
                 Pi ,
                  0.0000001,
                 10000);
    arma::mat tmpA = EM["A"];
    A = tmpA;
    arma::mat tmpPi = EM["Pi"];
   Pi = tmpPi;
     b_f = for_back(m, A, f0x, f1x, Pi);
    arma::mat gamma = b_f["gamma"];
    double s_poid = sum(gamma.col(1));
    for(int j = 0; j < m; j++){
     arma::vec  tmp_f1j = vec(m);
      for (int k = 0; k < m; k++){
        tmp_f1j(k) = exp(-0.5 * pow((x(k) - x(j) )/ h,2))  * gamma(k,1) ;
      }
      f1x(j) = sum(tmp_f1j)/ ( sqrt( 2 * M_PI) * s_poid * h);
    }
    diff(0) = abs(A - A_old).max();
    diff(1) = abs(Pi - Pi_old).max();
    diff(2) = abs(f1x - f1x_old).max();
    i++;
  } // End of while
  return Rcpp::List::create(Rcpp::Named("A") = A,
                            Rcpp::Named("Pi") = Pi,
                            Rcpp::Named("fw_bc_EM") = b_f,
                            Rcpp::Named("f1x") = f1x,
                            Rcpp::Named("i") = i);
}


//' Use EM algorithm to estimate the parameters A  Pi and f1 of an HMM
//'
//' @param m the number of positions (hypothesis)
//' @param alpha a matrix m * 2  containing the forward variables
//' @param beta a matrix m * 2  containing the backward variables
//' @param A a matrix 2 * 2 the transition probabilities
//' @param f0x a vector of the values of the density under the null hypothesis on the observations
//' @param f1x a vector of the values of the density under the alternative hypothesis on the observations
//' @param Pi a vector of the initial state probabilities
//' @param eps the value ta reach for the convergence
//' @param maxit integer, the maximum number of iteration
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
// [[Rcpp::export]]
Rcpp::List Em_tot(int m,
                 arma::mat A,
                 arma::vec Pi,
                 arma::vec f0x,
                 arma::vec f1x,
                 arma::vec x,
                 double eps,
                 int maxit,
                 double h) {
  arma::vec diff = vec(3);
  diff(0) = eps + 1;
  int i = 0;
  Rcpp::List b_f;
  while(((max(diff) > eps) & (i < maxit)))
  {
    arma::mat A_old = A;
    arma::vec Pi_old = Pi;
    arma::vec f1x_old = f1x;
    b_f = for_back(m, A, f0x, f1x, Pi);
    arma::mat ksi = b_f["ksi"];
    arma::mat gamma = b_f["gamma"];
    arma::rowvec sum_gamma = sum(gamma.rows(0,m-2));
    arma::rowvec sum_ksi = sum(ksi.rows(0,m-2));
    Pi = trans(gamma.row(1));
    for(int j = 0; j < 2; j++){
      for(int k = 0; k < 2; k++){
        A(j, k) = sum_ksi(j + 2 * k) / sum_gamma(j);
      }
    }
    double s_poid = sum(gamma.col(1));
    for(int j = 0; j < m; j++){
      arma::vec  tmp_f1j = vec(m);
      for (int k = 0; k < m; k++){
        tmp_f1j(k) = exp(-0.5 * pow((x(k) - x(j) )/ h,2))  * gamma(k,1) ;
      }
      f1x(j) = sum(tmp_f1j)/ ( sqrt( 2 * M_PI) * s_poid * h);
    }
    diff(0) = abs(A - A_old).max();
    diff(1) = abs(Pi - Pi_old).max();
    diff(2) = abs(f1x - f1x_old).max();
    i++;
  } // End of while
  
  b_f = for_back(m, A, f0x, f1x, Pi);
  return Rcpp::List::create(Rcpp::Named("A") = A,
                            Rcpp::Named("Pi") = Pi,
                            Rcpp::Named("fw_bc_EM") = b_f,
                            Rcpp::Named("f1x") = f1x,
                            Rcpp::Named("i") = i);
}



//' Use EM algorithm to estimate the parameters A and Pi of an HMM
//'
//' @param m the number of positions (hypothesis)
//' @param alpha a matrix m * 2  containing the forward variables
//' @param beta a matrix m * 2  containing the backward variables
//' @param A a matrix 2 * 2 the transition probabilities
//' @param f0x a vector of the values of the density under the null hypothesis on the observations
//' @param f1x a vector of the values of the density under the alternative hypothesis on the observations
//' @param Pi a vector of the initial state probabilities
//' @param eps the value ta reach for the convergence
//' @param maxit integer, the maximum number of iteration
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
// [[Rcpp::export]]
Rcpp::List Em_tot_01(int m,
                  arma::mat A,
                  arma::vec Pi,
                  arma::vec f0x,
                  arma::vec f1x,
                  arma::vec x,
                  double eps,
                  int maxit,
                  double h) {
  arma::vec diff = vec(4);
  diff(0) = eps + 1;
  int i = 0;
  Rcpp::List b_f;
  while(((max(diff) > eps) & (i < maxit)))
  {
    arma::mat A_old = A;
    arma::vec Pi_old = Pi;
    arma::vec f0x_old = f0x;
    arma::vec f1x_old = f1x;
    b_f = for_back(m, A, f0x, f1x, Pi);
    arma::mat ksi = b_f["ksi"];
    arma::mat gamma = b_f["gamma"];
    arma::rowvec sum_gamma = sum(gamma.rows(0,m-2));
    arma::rowvec sum_ksi = sum(ksi.rows(0,m-2));
    Pi = trans(gamma.row(1));
    for(int j = 0; j < 2; j++){
      for(int k = 0; k < 2; k++){
        A(j, k) = sum_ksi(j + 2 * k) / sum_gamma(j);
      }
    }
    double s_poid1 = sum(gamma.col(1));
    double s_poid0 = sum(gamma.col(0));
    for(int j = 0; j < m; j++){
      arma::vec  tmp_f1j = vec(m);
      arma::vec  tmp_f0j = vec(m);
      for (int k = 0; k < m; k++){
        tmp_f1j(k) = exp(-0.5 * pow((x(k) - x(j) )/ h,2))  * gamma(k,1) ;
        tmp_f0j(k) = exp(-0.5 * pow((x(k) - x(j) )/ h,2))  * gamma(k,0) ;
      }
      f1x(j) = sum(tmp_f1j)/ ( sqrt( 2 * M_PI) * s_poid1 * h);
      f0x(j) = sum(tmp_f0j)/ ( sqrt( 2 * M_PI) * s_poid0 * h);
    }
    diff(0) = abs(A - A_old).max();
    diff(1) = abs(Pi - Pi_old).max();
    diff(2) = abs(f1x - f1x_old).max();
    diff(3) = abs(f0x - f0x_old).max();
    i++;
  } // End of while
  
  b_f = for_back(m, A, f0x, f1x, Pi);
  return Rcpp::List::create(Rcpp::Named("A") = A,
                            Rcpp::Named("Pi") = Pi,
                            Rcpp::Named("fw_bc_EM") = b_f,
                            Rcpp::Named("f1x") = f1x,
                            Rcpp::Named("f0x") = f0x,
                            Rcpp::Named("i") = i);
}