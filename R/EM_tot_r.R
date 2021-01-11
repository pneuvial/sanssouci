#' Title
#'
#' @param m 
#' @param A 
#' @param Pi 
#' @param f0x 
#' @param f1x 
#' @param x 
#' @param eps 
#' @param maxit 
#' @param h 
#'
#' @return
#' @export
#'
#' @examples
#' m <- 200
#' A <- matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T)
#' Pi <- c(0.8,0.2)
#' theta <- sim_markov(m, Pi , A)
#' x <- rep(0, m)
#' x[theta == 0] <- rnorm(sum(theta ==0))
#' x[theta == 1] <- rnorm(sum(theta ==1), 2, 1)
#' f0x <- dnorm(x)
#' f1x <- dnorm(x, 2,1)
#' Em_tot_01_approx(m, A, Pi, f0x, f1x,
#'  x, eps= 0.0001, maxit=1000, h= 0.3)
Em_tot_01_approx<- function(m, A, Pi, f0x, f1x,
   x, eps, maxit, h) {
 diff <- rep(eps +1, 4)
   i <- 0
 
  browser()
  while(((max(diff) > eps) & (i < maxit)))
  {
    A_old = A
    Pi_old = Pi
    f0x_old = f0x
    f1x_old = f1x
    b_f = for_back(m, A, f0x, f1x, Pi)
    sum_gamma = colSums(b_f$gamma[-m,])
    sum_ksi = colSums(b_f$ksi[-m,])
    Pi = b_f$gamma[1,]
    A <- matrix(sum_ksi / rep(sum_gamma,2), byrow = FALSE, ncol = 2)
    d0 <- density(x, weights = b_f$gamma[,1]/sum(b_f$gamma[,1]))
    d1 <- density(x, weights = b_f$gamma[,2]/sum(b_f$gamma[,2]))
    f1x <- approx(d1$x,d1$y,x)$y
    f0x <- approx(d0$x,d0$y,x)$y
    
    diff[1] = max(abs(A - A_old))
    diff[2] = max(abs(Pi - Pi_old))
    diff[3] = max(abs(f1x - f1x_old))
    diff[4] = max(abs(f0x - f0x_old))
    
   
  }
  
  b_f = for_back(m, A, f0x, f1x, Pi)
  return (list(A= A,Pi = Pi,
                            fw_bc_EM= b_f,
                            f1x = f1x,
                           f0x = f0x,
                            i = i))
}
