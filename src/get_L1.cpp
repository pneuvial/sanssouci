#include "RcppArmadillo.h"

using namespace arma;

//' Simple matrix multiplication
//'
//' @param m the number of positions (hypothesis)
//' @param alpha a matrix m * 2  containing the forward variables
//' @param beta a matrix m * 2  containing the backward variables
//' @param A a matrix 2 * 2 the transition probabilities
//' @param f0x a vector of the values of the density under the null hypothesis on the observations
//' @param f1x a vector of the values of the density under the alternative hypothesis on the observations
//' @param i the position (hypothesis) for wich we want the posterior transition matrix.
//' @return Product of matrices
//' @export
//'
//' @examples
//' A <- matrix(1:9, 3, 3)
//' B <- matrix(11:19, 3, 3)
//' matrix_mult_cpp(A, B)
// [[Rcpp::export]]
arma::mat get_A_old(int m, arma::mat alpha,
          arma::mat beta,
          arma::mat A,
          arma::vec f0x,
          arma::vec f1x,
          int i) {
  double ai01;
  double ai11;
  double ai10;
  double ai00;
  double li0;
  double li1;
  double sum;
  arma::mat Ai(2,2);
  i = i -1;
  li0 = alpha(i - 1 , 0) * beta(i - 1 , 0) / (alpha(i - 1 , 0) * beta(i - 1 , 0) + alpha(i - 1 , 1) * beta(i - 1 , 1));
  li1 = alpha(i - 1 , 1) * beta(i - 1 , 1) / (alpha(i - 1 , 0) * beta(i - 1 , 0) + alpha(i - 1 , 1) * beta(i - 1 , 1));
  ai01 =  beta(i, 1) * alpha(i - 1, 0) * f1x[i] * A(0, 1);
  ai00 =  beta(i, 0) * alpha(i - 1, 0) * f0x[i] * A(0, 0);
  ai10 =  beta(i, 0) * alpha(i - 1, 1) * f0x[i] * A(1, 0);
  ai11 =  beta(i, 1) * alpha(i - 1, 1) * f1x[i] * A(1, 1);
  sum = ai00 + ai01 + ai10 + ai11;
  Ai(0,0) = ai00 / (li0*sum);
  Ai(0,1) = ai01 / (li0*sum);
  Ai(1,0) = ai10 / (li1*sum);
  Ai(1,1) = ai11 / (li1*sum);

  return Ai;
}


//' Simple matrix multiplication
//'
//' @param m the number of positions (hypothesis)
//' @param alpha a matrix m * 2  containing the forward variables
//' @param beta a matrix m * 2  containing the backward variables
//' @param A a matrix 2 * 2 the transition probabilities
//' @param f0x a vector of the values of the density under the null hypothesis on the observations
//' @param f1x a vector of the values of the density under the alternative hypothesis on the observations
//' @param i the position (hypothesis) for wich we want the posterior transition matrix.
//' @return Product of matrices
//' @export
//'
//' @examples
//' m <-  100
//' A <- matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T)
//'   rdata <- sim_hmm_2states(m, Pi, A, f0 = c(0,1), f1= c(2,1))
//'  x <- rdata$x
//'  f0x <- dnorm(x)
//'  f1x <- dnorm(x, 1,2)
//'   mod <- for_back(m, A, f0x, f1x, Pi)
//'     Pis_est <- lapply(2:m, function(i){
//'  get_A( m,alpha = mod$alpha, beta = mod$beta, A, f0x, 
//'         f1x, i = i)
//'  })
//' Pis_est2 <- lapply(2:m, function(i){
//'  get_A_old( m,alpha = mod$alpha, beta = mod$beta, A, f0x, 
//'         f1x, i = i)
//'  })
// [[Rcpp::export]]
arma::mat get_A(int m, arma::mat alpha,
                    arma::mat beta,
                    arma::mat A,
                    arma::vec f0x,
                    arma::vec f1x,
                    int i) {
  double ai01;
  double ai11;
  double ai10;
  double ai00;
  double li0;
  double li1;
  double sum;
  arma::mat Ai(2,2);
  i = i -1;
  li0 =  beta(i - 1 , 0) / (alpha(i - 1 , 0) * beta(i - 1 , 0) + alpha(i - 1 , 1) * beta(i - 1 , 1));
  li1 =  beta(i - 1 , 1) / (alpha(i - 1 , 0) * beta(i - 1 , 0) + alpha(i - 1 , 1) * beta(i - 1 , 1));
  ai01 =  beta(i, 1)  * f1x[i] * A(0, 1);
  ai00 =  beta(i, 0) * f0x[i] * A(0, 0);
  ai10 =  beta(i, 0)  * f0x[i] * A(1, 0);
  ai11 =  beta(i, 1) * f1x[i] * A(1, 1);
  sum = ai00 * alpha(i - 1, 0) + ai01 * alpha(i - 1, 0)+ ai10  * alpha(i - 1, 1)+ ai11 * alpha(i - 1, 1);
  Ai(0,0) = ai00 / (li0*sum);
  Ai(0,1) = ai01 / (li0*sum);
  Ai(1,0) = ai10 / (li1*sum);
  Ai(1,1) = ai11 / (li1*sum);
  
  return Ai;
}



//' Simple matrix multiplication
//'
//' @param m the number of positions (hypothesis)
//' @param alpha a matrix m * 2  containing the forward variables
//' @param beta a matrix m * 2  containing the backward variables
//' @param A a matrix 2 * 2 the transition probabilities
//' @param f0x a vector of the values of the density under the null hypothesis on the observations
//' @param f1x a vector of the values of the density under the alternative hypothesis on the observations
//'
//' @return Product of matrices
//' @export
//'
//' @examples
//' A <- matrix(1:9, 3, 3)
//' B <- matrix(11:19, 3, 3)
//' matrix_mult_cpp(A, B)
// [[Rcpp::export]]
arma::sp_mat get_L1(arma::mat A, int m,
              arma::mat alpha,
              arma::mat beta,
              arma::vec f0x,
              arma::vec f1x) {
  arma::sp_mat full(m,m);
  arma::sp_mat A_bdiag(2*m, 2*m);
  arma::mat A1;
  A1 = get_A( m,  alpha,
              beta,
              A,
              f0x,
              f1x, m);
  full(m-1, m-2) = A1(0,0);
  A_bdiag(span( 2*(m-1), 2*m-1), span( 2*(m-1), 2*m-1)) = A1;
  for(int i = m-1; i > 1 ; --i) {
    int num;
    num = 2+ 2*(m-i);
    arma::sp_mat old_A_bdiag(num,num);

    arma::mat Ai(2,2);
    Ai = get_A( m,  alpha,
                beta,
                A,
                f0x,
                f1x, i );
    old_A_bdiag(span(0,1), span(0,1)) = speye(2,2);
    old_A_bdiag(span(2, num-1), span(2, num-1)) = A_bdiag(span( 2*i, 2*m-1), span( 2*i, 2*m-1));
    arma::sp_mat A_mom ( num, num);
    for ( int j = 0; j <  m - i +1 ; ++j){
      A_mom(span(2*j, 2*j+1),span(2*j, 2*j+1)) = Ai;
    }
    A_bdiag(span( 2*(i-1), 2*m-1), span( 2*(i-1), 2*m-1)) = A_mom * old_A_bdiag;
    arma::sp_mat pti(2*(m -i)+2,2*(m -i)+2);
    pti = A_bdiag(span( 2*(i-1), 2*m-1), span( 2*(i-1), 2*m-1));
    arma::uvec a = regspace<uvec>(0, 2, 2*(m -i));
    arma::vec  v(pti.diag());
    full(span(i-1, m-1),i-2) =v.elem(a);
  }
  full.diag()  += 1;
  return(full);
}

//' New way of finding Bin ! (Now A)
//'
//' @param A a matrix 2 * 2 the transition probabilities
//' @param m the number of positions (hypothesis)
//' @param alpha a matrix m * 2  containing the forward variables
//' @param beta a matrix m * 2  containing the backward variables
//' @param f0x a vector of the values of the density under the null hypothesis on the observations
//' @param f1x a vector of the values of the density under the alternative hypothesis on the observations
//'
//' @return Product of matrices
//' @export
//'
//' @examples
//' A <- matrix(1:9, 3, 3)
//' B <- matrix(11:19, 3, 3)
//' matrix_mult_cpp(A, B)
// [[Rcpp::export]]
Rcpp::List getA01( int m,
                   arma::vec li0,
                   arma::vec f0x,
                   arma::vec f1x,
                   Rcpp::List Pis){

   arma::sp_mat B0(m, m + 1);
   arma::sp_mat B1(m, m + 1);
   B0(0, span(1, m )) += li0(0);
   B1(0, span(0, m )) += 1 - li0(0);
   B0(span(0, m - 1), 0) += 0;
   for ( int i = 1; i <  m  ; ++i){
      arma::mat Pii = Pis(i - 1);
      B1(i, 0) = B1(i - 1, 0) * Pii(1, 1);
      for ( int j = 1; j <  m  + 1  ; ++j){
         B0(i, j) = B0(i  -1, j -1) * Pii(0,0) + B1(i  -1, j -1) * Pii(1,0);
         B1(i, j) = B0(i  -1, j) * Pii(0,1) + B1(i  -1, j ) * Pii(1,1);
      }
   }
   return Rcpp::List::create(Rcpp::Named("A1")=B1,
                             Rcpp::Named("A0")= B0);
}



//' New way of finding Bin ! (Now A)
//'
//' @param A a matrix 2 * 2 the transition probabilities
//' @param m the number of positions (hypothesis)
//' @param alpha a matrix m * 2  containing the forward variables
//' @param beta a matrix m * 2  containing the backward variables
//' @param f0x a vector of the values of the density under the null hypothesis on the observations
//' @param f1x a vector of the values of the density under the alternative hypothesis on the observations
//'
//' @return Product of matrices
//' @export
//'
//' @examples
//' A <- matrix(1:9, 3, 3)
//' B <- matrix(11:19, 3, 3)
//' matrix_mult_cpp(A, B)
// [[Rcpp::export]]
double getbound( int m, double alpha, 
                   arma::vec li0,
                   arma::vec f0x,
                   arma::vec f1x,
                   Rcpp::List Pis){
  
  arma::sp_mat B0(m, m + 1);
  arma::sp_mat B1(m, m + 1);
  B0(0, span(1, m )) += li0(0);
  B1(0, span(0, m )) += 1 - li0(0);
  B0(span(0, m - 1), 0) += 0;
  for ( int i = 1; i <  m  ; ++i){
    arma::mat Pii = Pis(i - 1);
    B1(i, 0) = B1(i - 1, 0) * Pii(1, 1);
  }
 // for ( int j = 1; j <  m  + 1  ; ++j){

  double risk = B0(m-1, 0) + B1(m-1, 0); 
  int j =1;
     while(risk < (1 -alpha)){
       
    for ( int i = 1; i <  m  ; ++i){
      arma::mat Pii = Pis(i - 1);
      B0(i, j) = B0(i  -1, j -1) * Pii(0,0) + B1(i  -1, j -1) * Pii(1,0);
      B1(i, j) = B0(i  -1, j) * Pii(0,1) + B1(i  -1, j ) * Pii(1,1);
    }
    risk =  B0(m-1, j) + B1(m-1, j); 
    ++j; 
  }
  return((j-1) );
}


//' New way of finding Bin ! (Now A)
//'
//' @param A a matrix 2 * 2 the transition probabilities
//' @param m the number of positions (hypothesis)
//' @param alpha a matrix m * 2  containing the forward variables
//' @param beta a matrix m * 2  containing the backward variables
//' @param f0x a vector of the values of the density under the null hypothesis on the observations
//' @param f1x a vector of the values of the density under the alternative hypothesis on the observations
//'
//' @return Product of matrices
//' @export
//'
//' @examples
//' m <-  100
//' A <- matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T)
//'   rdata <- sim_hmm_2states(m, Pi, A, f0 = c(0,1), f1= c(2,1))
//'  x <- rdata$x
//'  f0x <- dnorm(x)
//'  f1x <- dnorm(x, 1,2)
//'   mod <- for_back(m, A, f0x, f1x, Pi)
//'     Pis_est <- lapply(2:m, function(i){
//'  get_A( m,alpha = mod$alpha, beta = mod$beta, A, f0x, 
//'         f1x, i = i)
//'  })
//'  alpha <- 0.1
//' getIC(length(x), alpha, mod$gamma[,1], f0x, f1x, Pis_est)
//' quant <- get_quantiles(sel = 1:m, li0 = mod$gamma[, 1], 
//'  Pis = Pis_est, f0x = f0x, f1x = f1x)
//'  borne(type_borne = "HMM", sel = 1:m, a = quant, alpha )
//'  borne(type_borne = "HMM_small", sel = 1:m, a = quant, alpha)
// [[Rcpp::export]]
arma::vec getIC( int m, double alpha, 
                 arma::vec li0,
                 arma::vec f0x,
                 arma::vec f1x,
                 Rcpp::List Pis){
  
  arma::sp_mat B0(m, m + 1);
  arma::sp_mat B1(m, m + 1);
  B0(0, span(1, m )) += li0(0);
  B1(0, span(0, m )) += 1 - li0(0);
  B0(span(0, m - 1), 0) += 0;
  for ( int i = 1; i <  m  ; ++i){
    arma::mat Pii = Pis(i - 1);
    B1(i, 0) = B1(i - 1, 0) * Pii(1, 1);
  }
  // for ( int j = 1; j <  m  + 1  ; ++j){
  
  double risk = B0(m-1, 0) + B1(m-1, 0); 
  int j =1;
  arma::vec bounds(2);
  while(risk < alpha){
    
    for ( int i = 1; i <  m  ; ++i){
      arma::mat Pii = Pis(i - 1);
      B0(i, j) = B0(i  -1, j -1) * Pii(0,0) + B1(i  -1, j -1) * Pii(1,0);
      B1(i, j) = B0(i  -1, j) * Pii(0,1) + B1(i  -1, j ) * Pii(1,1);
    }
    risk =  B0(m-1, j) + B1(m-1, j); 
    ++j; 
  }
  bounds(0) = j -2;
  
  while(risk < (1 -alpha)){
    for ( int i = 1; i <  m  ; ++i){
    arma::mat Pii = Pis(i - 1);
    B0(i, j) = B0(i  -1, j -1) * Pii(0,0) + B1(i  -1, j -1) * Pii(1,0);
    B1(i, j) = B0(i  -1, j) * Pii(0,1) + B1(i  -1, j ) * Pii(1,1);
  }
  risk =  B0(m-1, j) + B1(m-1, j); 
  ++j; 
}
  bounds(1) = j -1;
  return(bounds);
}

//' New way of finding Bin ! (Now A)
//'
//' @param A a matrix 2 * 2 the transition probabilities
//' @param m the number of positions (hypothesis)
//' @param alpha a matrix m * 2  containing the forward variables
//' @param beta a matrix m * 2  containing the backward variables
//' @param f0x a vector of the values of the density under the null hypothesis on the observations
//' @param f1x a vector of the values of the density under the alternative hypothesis on the observations
//'
//' @return Product of matrices
//' @export
//'
//' @examples
//' m <-  100
//' A <- matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T)
//'   rdata <- sim_hmm_2states(m, Pi, A, f0 = c(0,1), f1= c(2,1))
//'  x <- rdata$x
//'  f0x <- dnorm(x)
//'  f1x <- dnorm(x, 1,2)
//'   mod <- for_back(m, A, f0x, f1x, Pi)
//'     Pis_est <- lapply(2:m, function(i){
//'  get_A( m,alpha = mod$alpha, beta = mod$beta, A, f0x, 
//'         f1x, i = i)
//'  })
//'  alpha <- 0.1
//' getIC(length(x), alpha, mod$gamma[,1], f0x, f1x, Pis_est)
//' quant <- get_quantiles(sel = 1:m, li0 = mod$gamma[, 1], 
//'  Pis = Pis_est, f0x = f0x, f1x = f1x)
//'  borne(type_borne = "HMM", sel = 1:m, a = quant, alpha )
//'  borne(type_borne = "HMM_small", sel = 1:m, a = quant, alpha)
// [[Rcpp::export]]
arma::vec getprobs( int m, arma::vec prob, 
                 int size_prob,
                 arma::vec petit_grand,
                 arma::vec li0,
                 arma::vec f0x,
                 arma::vec f1x,
                 Rcpp::List Pis){
  
  arma::sp_mat B0(m, m + 1);
  arma::sp_mat B1(m, m + 1);
  B0(0, span(1, m )) += li0(0);
  B1(0, span(0, m )) += 1 - li0(0);
  B0(span(0, m - 1), 0) += 0;
  for ( int i = 1; i <  m  ; ++i){
    arma::mat Pii = Pis(i - 1);
    B1(i, 0) = B1(i - 1, 0) * Pii(1, 1);
  }
  // for ( int j = 1; j <  m  + 1  ; ++j){
  
  double risk = B0(m-1, 0) + B1(m-1, 0); 
  int j =1;
  arma::vec bounds(size_prob);
  for ( int p = 0; p < size_prob  ; ++p){
  while(risk < prob(p)){
    
    for ( int i = 1; i <  m  ; ++i){
      arma::mat Pii = Pis(i - 1);
      B0(i, j) = B0(i  -1, j -1) * Pii(0,0) + B1(i  -1, j -1) * Pii(1,0);
      B1(i, j) = B0(i  -1, j) * Pii(0,1) + B1(i  -1, j ) * Pii(1,1);
    }
    risk =  B0(m-1, j) + B1(m-1, j); 
    ++j; 
  }
  bounds(p) = j -petit_grand(p) -1;
  }
  return(bounds);
}

