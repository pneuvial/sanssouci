#' Title
#'
#' @param x_obs
#' @param f0x
#' @param x
#' @param h
#'
#' @return
#' @export
#'
#' @examples
f_hat <- function(x_obs, x, h = 0.3){
  sum(abs(x_obs - x) < h) / (2 * length(x_obs) * h)
}

#' Title
#'
#' @param f0x
#' @param f_hatx
#' @param pi0_hat
#'
#' @return
#' @export
#'
#' @examples
f1x_hat <- function(f0x, f_hatx, pi0_hat){
 a <- (f_hatx - pi0_hat*f0x) /(1 -pi0_hat)
 a[a<0]<-0
 return(a)
 }
