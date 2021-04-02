#' Title
#'
#' @param sel 
#' @param li0 
#' @param Pis 
#' @param f0x 
#' @param f1x 
#'
#' @return
#' @export
#'
#' @examples
get_quantiles <- function(sel, li0, Pis, f0x, f1x){
  if(length(sel)< 2) {
    return(NA)
  }else{
    Pis_sel <- lapply(1:(length(sel)-1), function(i){
      Reduce("%*%", Pis[sel[i]: (sel[i +1]-1)])
    })
    l_sel <- length(sel)
    A01 <- getA01(m = l_sel, 
                  li0 = li0[sel], 
                  f0x[sel],
                  f1x[sel],
                  Pis_sel)  
    a <- A01$A1[
      l_sel, ] + A01$A0[l_sel, ]
    return(a)
  }
}


#' Title
#'
#' @param sel must be ordered!
#' @param li0 
#' @param Pis 
#' @param f0x 
#' @param f1x 
#'
#' @return
#' @export
#'
#' @examples
#'  m <-  100
#'  Pi <- c(0.8,0.2)
#'  A <- matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T)
#'    rdata <- sim_hmm_2states(m, Pi, A, f0 = c(0,1), f1= c(2,1))
#'   x <- rdata$x
#'   f0x <- dnorm(x)
#'   f1x <- dnorm(x, 1,2)
#'    mod <- for_back(m, A, f0x, f1x, Pi)
#'      Pis_est <- lapply(2:m, function(i){
#'   get_A( m,alpha = mod$alpha, beta = mod$beta, A, f0x, 
            #'          f1x, i = i)
#'   })
#'   alpha <- 0.1
#'   sel <- sample(1:m, m/2)
#'   sel <- sel[order(sel)]
#'  get_IC(sel = sel, li0 = mod$gamma[, 1], 
#'   Pis = Pis_est, f0x = f0x, f1x = f1x, alpha)
#'  quant <- get_quantiles(sel = sel, li0 = mod$gamma[, 1], 
#'   Pis = Pis_est, f0x = f0x, f1x = f1x)
#'   borne(type_borne = "HMM", sel = sel, a = quant, alpha )
#'   borne(type_borne = "HMM_small", sel = sel, a = quant, alpha)
get_IC <- function(sel, li0, Pis, f0x, f1x, alpha){
  if(length(sel)< 2) {
    return(NA)
  }else{
    Pis_sel <- lapply(1:(length(sel)-1), function(i){
      Reduce("%*%", Pis[sel[i]: (sel[i +1]-1)])
    })
    l_sel <- length(sel)
    IC <- getIC(m = l_sel, 
                alpha,
                  li0 = li0[sel], 
                  f0x[sel],
                  f1x[sel],
                  Pis_sel)  
    
    return(IC)
  }
}


#' Title
#'
#' @param sel 
#' @param li0 
#' @param Pis 
#' @param f0x 
#' @param f1x 
#' @param probs 
#'
#' @return
#' @export
#'
#' @examples
get_probs <- function(sel, li0, Pis, f0x, f1x, probs){
  if(length(sel)< 2) {
    return(NA)
  }else{
    Pis_sel <- lapply(1:(length(sel)-1), function(i){
      Reduce("%*%", Pis[sel[i]: (sel[i +1]-1)])
    })
    l_sel <- length(sel)
    probs_ord <- sort(probs)
    petit_grand <- rep(0, length(probs))
    petit_grand[probs_ord <0.5] <- 1
    Probs <- getprobs( length(sel),  
                       probs_ord, 
                       length(probs),
                       petit_grand,
                       li0 = li0[sel], 
                       f0x[sel],
                       f1x[sel],
                       Pis_sel)
    Probs_ord <- Probs[order(order(probs))]
    return(Probs_ord)
  }
}
