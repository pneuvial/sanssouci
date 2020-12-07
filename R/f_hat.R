
#' Non parametric density estimation
#'
#' @param x_obs numeric vector : the whole observed values
#' @param x numeric, the value where the density is estimated
#' @param h numeric, the size of the window
#' @param K function, the kernel to use 
#'
#' @return
#' @export
#'
#' @examples
f_hatK <- function(x_obs, x, h = 0.3, K){
  sum(K((x_obs - x) / h)) / ( length(x_obs) * h)
}


#' Naive f1 estimate
#'
#' @param f0x numeric vector.  The value of the density for the different obersvation under the null hypothesis.
#' @param f_hatx numeric vector.  The value of the estimated density of the different obersvation.
#' @param pi0_hat numeric, the estimated proportion of the null hypothesis. 
#'
#' @return numeric vector.  The value of the estimated density for the different obersvation under the alternatives.
#' @export 
#'
#' @examples
f1x_hat <- function(f0x, f_hatx, pi0_hat){
 a <- (f_hatx - pi0_hat*f0x) /(1 -pi0_hat)
 a[a<0]<-0
 return(a)
 }
