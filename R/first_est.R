#' Title
#'
#' @param x 
#' @param h 
#' @param m0_init 
#' @param sd0_init 
#' @param df_init 
#' @param norm_init 
#' @param max_pi0 
#' @param f0_known 
#' @param f0x_est 
#' @param pval 
#' @param plot 
#' @param size_plot 
#'
#' @return
#' @export
#'
#' @examples
#' m <- 2000
#' theta <- sim_markov(m, Pi = c(0.8,0.2), A = matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T))
#' x <- rep(0, m)
#'  x[theta == 0] <- rnorm(sum(theta ==0))
#'  x[theta == 1] <- rnorm(sum(theta ==1), 2, 1)
#'  Estimation(x, m0_init = 0, sd0_init = 1, norm_init = TRUE, plot= TRUE) 
Estimation <- function(x,  h =0.3,
                     m0_init, sd0_init, df_init, norm_init, max_pi0= 0.99999, 
                     f0_known = TRUE, f0x_est = NULL, pval = NULL, 
                     plot = FALSE, size_plot= min(10000, length(x)), 
                     approx = TRUE){
  m <- length(x)
  
if(is.null(f0x_est)){
    if(norm_init){
    f0x_est <- dnorm(x, m0_init, sd0_init)
  }else{
    f0x_est <- dt(x, df_init)
  }

}
  if(is.null(pval)){
    
    if(norm_init){
      pval <- 2 * (1 - pnorm(abs(x), m0_init, sd0_init))
    }else{
      pval <- 2 * (1 - pt(abs(x), df_init)) 
    }
    
  }
pi0_hat <- max(min(sum(pval > 0.8) / (m * 0.2), max_pi0), 0.6)
if(approx){
  d <- density(x,bw = h)
  f_hatx <- approx(d$x,d$y,x)$y
}else {
  f_hatx <- x %>%
    map_dbl( ~f_hatK(x, ., h = h,K) )
}

f1x_est <-  f1x_hat(f0x_est, f_hatx, pi0_hat)
f1x_est[f1x_est <= 0] <- min(f0x_est)
mini <- max(0.6, ((1 + max_pi0) * pi0_hat -max_pi0) / pi0_hat)
a <-runif(1, mini, max_pi0)
b <- 1 - a
c <- pi0_hat * b / (1 - pi0_hat)
d <- 1 - c
Em <- Em(m, A = matrix(c(a, b, c, d), byrow = TRUE, ncol=2),
                Pi= c(pi0_hat, 1 -pi0_hat),  f0x = f0x_est, f1x = f1x_est,
                x, eps = 0.0001,
                maxit =1000, h = h, f0_known, approx = approx)

p <- NULL
if(plot){
  nb_plot <- round(seq(1, length(x), length.out = size_plot))
  
if(!f0_known){
    p<- tibble(x =x[nb_plot],f0x_init = f0x_est[nb_plot], f0x_est = Em$f0x[nb_plot], 
           f1x_init = f1x_est[nb_plot], f1x_est = Em$f1x[nb_plot] ) %>% 
      gather(-x, key = "key", value = "value") %>% 
      ggplot(aes(x = x, y = value, color = key)) + geom_line() + theme_bw()
  }else{
    p<- tibble(x =x[nb_plot],f0x_init = f0x_est[nb_plot],  
               f1x_init = f1x_est[nb_plot], f1x_est = Em$f1x[nb_plot] ) %>% 
      gather(-x, key = "key", value = "value") %>% 
      ggplot(aes(x = x, y = value, color = key)) + geom_line() + theme_bw()
    
  }
  print(p)
}

return(list(Em =Em, plot = p))

}
