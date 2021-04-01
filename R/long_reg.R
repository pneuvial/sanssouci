#' Title
#'
#' @param viterbi 
#' @param min_size 
#'
#' @return
#' @export
#'
#' @examples
long_reg <- function(viterbi, min_size){
  nb_0 <- cumsum( 1 - viterbi)
  Rk_s <- split(which(viterbi == 1), nb_0[viterbi == 1])
  size <- sapply(Rk_s, length)
  Rk_s <- Rk_s[which(size >= min_size)]
  unlist(Rk_s)
}



#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
K <- function(x){
  1 / sqrt(2 * pi) * exp(-x^2/2 )
}

#' Title
#'
#' @param x_obs 
#' @param x 
#' @param h 
#' @param K 
#'
#' @return
#' @export
#'
#' @examples
f_hatK <- function(x_obs, x, h = 0.3, K){
  sum(K((x_obs - x) / h)) / ( length(x_obs) * h)
}


#' Title
#'
#' @param A 
#' @param Pi 
#' @param x_from 
#' @param prob1 
#' @param h 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
sim_hmm_from_weightkde <- function( A, Pi,  x_from, prob1, h,n ){
  
  if(any(is.na(Pi))| max(Pi) > 1 | min(Pi)<0 | max(A) > 1 | min(A) < 0){
    stop('Pb de PIIIII')
  }
  m <- length(x_from)
  theta <- rep(0, m)
  x <- rep(0, m)
  
  ## generating theta :
  theta[1] <- rbinom(1, 1, Pi[2])
  for (i in 2:m)
  { theta[i] <- (1-theta[i-1])*rbinom(1, 1, A[1, 2]) + 
    theta[i-1]*rbinom(1, 1, A[2, 2])
  }
  
  
  for (ind in 0:1)
  { nb.ind <- sum(theta==ind)
  if (nb.ind>0) {
    if (ind == 0){
      x[theta==ind] <-  sample(x_from, nb.ind, replace = TRUE, prob  = 1- prob1) + 
        rnorm(nb.ind, sd = h)
    } else{
      x[theta==ind] <-  sample(x_from, nb.ind, replace = TRUE, prob  = prob1) +
        rnorm(nb.ind, sd = h)
    }
    
  }
  
  }
  
  data<-list(theta=theta, x=x)
  return (data)
  
}
