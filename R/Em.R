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
#' @param f0_known 
#' @param approx 
#'
#' @return
#' @export
#'
#' @examples
Em <- function(m, A ,
               Pi,  f0x, f1x,
               x, eps = 0.0001,
               maxit =1000, h = 0.3, f0_known, approx = TRUE){
  if(approx){
  if(f0_known){
    Em_tot_approx(m, A ,
              Pi,  f0x, f1x,
              x, eps,
              maxit, h)
  }else{
    Em_tot_01_approx(m, A ,
              Pi,  f0x, f1x,
              x, eps,
              maxit, h)   
  }
  }else{
    if(f0_known){
      Em_tot(m, A ,
             Pi,  f0x, f1x,
             x, eps,
             maxit, h)
    }else{
      Em_tot_01(m, A ,
                Pi,  f0x, f1x,
                x, eps,
                maxit, h)   
    }  
    }

}