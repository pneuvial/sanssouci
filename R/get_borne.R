#' Title
#'
#' @param B the matrix defined in the paper
#' @param alpha the levels 
#'
#' @return
#' @export
#'
#' @examples
#'  m     <-  100
#'  A     <- matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T)
#'  f0    <- c(0, 1)
#'  f1    <- c(2, 1)
#'  Pi    <- c( 0.9, 0.1)
#'  rdata <- sim_hmm_2states(m, Pi, A, f0, f1)
#'  x     <- rdata$x
#'  f0x   <- dnorm(x, 0, 1)
#'  f1x   <- dnorm(x, SNR, 1)
#'  mod   <- for_back(m, A, f0x, f1x, Pi)
#' viter_log <- viterbi_log(m,log(A),log(f0x),log(f1x), log(Pi))
#' Rk_s <- get_Rks(viterbi_log, min_size =3)
#' L <- get_L1(A, m,fw_bc_oracle$alpha, fw_bc_oracle$beta, f0x,  f1x)
#' Zeta <- get_zeta_ks(Rk_s, al = 0.2, for_back = fw_bc_oracle, L =L, f0x,  f1x)
#' get_phborne(S = 1:m, Rk_s, Zeta)
get_borne <- function(B, alpha){

  proba <- apply(B, 2, sum)
  if(min(proba)<= alpha){
 min(which(proba <= alpha))  - 1
  }else{
      ncol(B)
    }
}



#' Title
#'
#' @param viterbi the  path of 0/1 obtained using viterbi algorithm
#' @param min_size minimum size of a set "Rk"
#' @return
#' @export
#'
#' @examples
#'  m <-  100
#'  A <- matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T)
#'  f0 <- c(0, 1)
#'  f1 <- c(2, 1)
#'  Pi <- c( 0.9, 0.1)
#'  rdata <- simulate.data.hmm.2states(m, Pi, A, f0, f1)
#'  x <- rdata$x
#'  theta <- rdata$theta
#'  mod <- for_back(m, A, f0x, f1x, Pi)
#' viter_log <- viterbi_log(m,log(A),log(f0x),log(f1x), log(Pi))
#' Rk_s <- get_Rks(viterbi_log, min_size =3)
#' L <- get_L1(A, m,fw_bc_oracle$alpha, fw_bc_oracle$beta, f0x,  f1x)
#' Zeta <- get_zeta_ks(Rk_s, al = 0.2, for_back = fw_bc_oracle, L =L, f0x,  f1x)
#' get_phborne(S = 1:m, Rk_s, Zeta)
get_Rks <- function(viterbi, min_size = 1){
  nb_0 <- cumsum( 1 - viterbi)
 Rk_s <- split(which(viterbi == 1), nb_0[viterbi == 1])
 if(min_size !=1){
   size <- sapply(Rk_s, length)
   Rk_s <- Rk_s[which(size >= min_size)]
 }
 return(Rk_s)
}






#' Title
#'
#' @param R_ks list containing the position for the different Rk
#' @param al confidence level alpha
#' @param for_back result from the forward backward algorithm
#' @param L the L matrix (L[i,j] is the probability of theta_j = 0 knowing theta_i = 0 and X)
#' @param f0x the values of the density under the null hypothesis for the different points
#' @param f1x the values of the density under the alternative hypothesis for the different points
#'
#' @return
#' @export
#'
#' @examples
#'  m <-  100
#'  A <- matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T)
#'  f0 <- c(0, 1)
#'  f1 <- c(2, 1)
#'  Pi <- c( 0.9, 0.1)
#'  rdata <- simulate.data.hmm.2states(m, Pi, A, f0, f1)
#'  x <- rdata$x
#'  theta <- rdata$theta
#'  mod <- for_back(m, A, f0x, f1x, Pi)
#' viter_log <- viterbi_log(m,log(A),log(f0x),log(f1x), log(Pi))
#' Rk_s <- get_Rks(viterbi_log, min_size =3)
#' L <- get_L1(A, m,fw_bc_oracle$alpha, fw_bc_oracle$beta, f0x,  f1x)
#' Zeta <- get_zeta_ks(Rk_s, al = 0.2, for_back = fw_bc_oracle, L =L, f0x,  f1x)
#' get_phborne(S = 1:m, Rk_s, Zeta)
get_zeta_ks <- function(R_ks, al, for_back, Pis, f0x,  f1x){
  K <- length(R_ks)
  sapply(R_ks, function(R_k){
    if(length(R_k)>1){
      quant <- get_quantiles(sel = R_k, li0 =for_back$gamma[,1],
                             Pis = Pis, f0x = f0x, f1x = f1x)
      borne(type_borne = "HMM", sel = R_k, a = quant, alpha = al / K)
    }else{
      for_back$gamma[R_k,1]<= al/K
    }
  })
}


#' Title
#'
#' @param S
#' @param Rk_s
#' @param zeta_ks
#'
#' @return
#' @export
#'
#' @examples
get_phborne <- function(S, Rk_s, zeta_ks){
  K <- length(Rk_s)
  brns <- sapply(1:K, function(i){
    min(length(intersect(S, Rk_s[[i]])), zeta_ks[i])
  })
  sum(brns) + sum(!S %in% unlist(Rk_s))
}
