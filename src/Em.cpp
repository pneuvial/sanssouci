// #include "RcppArmadillo.h"
//
// using namespace arma;
//
// // This is a simple example of exporting a C++ function to R. You can
// // source this function into an R session using the Rcpp::sourceCpp
// // function (or via the Source button on the editor toolbar). Learn
// // more about Rcpp at:
// //
// //   http://www.rcpp.org/
// //   http://adv-r.had.co.nz/Rcpp.html
// //   http://gallery.rcpp.org/
// //
//
//
// //' Simple matrix multiplication
// //'
// //' @param A Matrix
// //' @param B Matrix
// //'
// //' @return Product of matrices
// //' @export
// //'
// //' @examples
// //' A <- matrix(1:9, 3, 3)
// //' B <- matrix(11:19, 3, 3)
// //' matrix_mult_cpp(A, B)
// // [[Rcpp::export]]
// Rcpp::List for_back(int m,
//                     arma::mat A,
//                     arma::vec f0x,
//                     arma::vec f1x,
//                     arma::vec Pi) {
//   arma::mat alpha(m, 2);
//   arma::vec c0(m);
//
//   c0(0) = 1 / (Pi(0)*f0x(0) + Pi(1)*f1x(0));
//   alpha(0, 0) = Pi(0)*f0x(0) * c0(0);
//   alpha(0, 1) = Pi(1)*f1x(0) * c0(0);
//
//
//   for (int i = 1; i < m ; ++i)
//   {
//
//     alpha(i, 0)  = f0x(i) * (alpha(i - 1, 0) * A(0, 0) + alpha(i - 1, 1) * A(1, 0));
//     alpha(i, 1)  = f1x(i) * (alpha(i - 1, 0) * A(0, 1) + alpha(i -1, 1) * A(1, 1));
//     c0(i) = 1 / ( alpha(i, 0) +  alpha(i, 1));
//     alpha(i, 0) = alpha(i, 0) * c0(i);
//     alpha(i, 1) = alpha(i, 1) * c0(i);
//   }
//   arma::mat beta(m, 2);
//   beta(m - 1, 0) = 1;
//   beta(m - 1, 1) = 1;
//
//   beta.row(m - 1) = c0(m - 1) * beta.row(m - 1);
//   arma::mat gamma(m, 2);
//   arma::mat ksi(m, 4);
//   for (int i = m - 2; i > -1; --i)
//   {
//     beta (i, 0) = c0(i) * (A(0, 0) * f0x(i + 1) * beta(i + 1, 0) + A(0, 1) * f1x(i + 1) * beta(i + 1, 1));
//     beta (i, 1) = c0(i) * (A(1, 0) * f0x(i + 1) * beta(i + 1, 0) + A(1, 1) * f1x(i + 1) * beta(i + 1, 1));
//     gamma(i, 0) = alpha(i, 0) * beta(i ,0) / ( alpha(i, 0) * beta(i ,0) + alpha(i, 1) * beta(i ,1));
//     gamma(i, 1) = alpha(i, 1) * beta(i ,1) / ( alpha(i, 0) * beta(i ,0) + alpha(i, 1) * beta(i ,1));
//     double denom =0;
//     for(int j = 0; j < 2; j++){
//       for(int k = 0; k < 2; k++){
//         double fx = (1-k)*f0x(i+1)+k*f1x(i+1);
//         denom =    alpha(i, j) * A(j, k) * fx * beta(i + 1, k) + denom;
//       }
//     }
//
//     for(int j = 0; j < 2; j++){
//       ksi(i, j) = alpha(i, j) * A(j, 0) * f0x(i + 1) * beta(i + 1, 0) /  denom;
//       ksi(i, j  + 2) = alpha(i, j) * A(j, 1) * f1x(i + 1) * beta(i + 1, 1) / denom;
//     }
//   }
//
//
//   return Rcpp::List::create(Rcpp::Named("alpha")=alpha,
//                             Rcpp::Named("beta")= beta,
//                             Rcpp::Named("gamma")= gamma,
//                             Rcpp::Named("ksi")= ksi);
// }
//
//
//
//
// //' Simple matrix multiplication
// //'
// //' @param A Matrix
// //' @param B Matrix
// //'
// //' @return Product of matrices
// //' @export
// //'
// //' @examples
// //' A <- matrix(1:9, 3, 3)
// //' B <- matrix(11:19, 3, 3)
// //' matrix_mult_cpp(A, B)
// // [[Rcpp::export]]
// Rcpp::List Em_hmm(int m,
//                   arma::mat A,
//                   arma::vec f0x,
//                   arma::vec f1x,
//                   arma::vec Pi,
//                   double eps,
//                   int maxit) {
//   arma::vec diff = vec(2);
//   diff(0) = eps + 1;
//   int i = 0;
//   while(max(diff) > eps & i < maxit )
//   {
//     arma::mat A_old = A;
//     arma::vec Pi_old = Pi;
//     Rcpp::List b_f = for_back(m, A, f0x, f1x, Pi);
//     arma::mat ksi = b_f["ksi"];
//     arma::mat gamma = b_f["gamma"];
//     arma::rowvec sum_gamma = sum(gamma.rows(0,m-2));
//     arma::rowvec sum_ksi = sum(ksi.rows(0,m-2));
//     Pi = trans(gamma.row(1));
//     for(int j = 0; j < 2; j++){
//       for(int k = 0; k < 2; k++){
//         A(j, k) = sum_ksi(j + 2 * k) / sum_gamma(j);
//       }
//     }
//     diff(0) = abs(A - A_old).max();
//     diff(1) = max(Pi - Pi_old);
//     i++;
//   }
//   return Rcpp::List::create(Rcpp::Named("A") = A,
//                             Rcpp::Named("Pi") = Pi,
//                             Rcpp::Named("i") = i);
// }
//
//
