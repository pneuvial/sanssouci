#' Title
#'
#' @param type_borne 
#' @param sel 
#' @param a 
#' @param alpha 
#' @param m 
#' @param pval 
#' @param C 
#' @param ZL 
#' @param leaf_list 
#'
#' @return
#' @export
#'
#' @examples
borne <- function( type_borne, sel, a, alpha, m, 
                   pval,
                   C, ZL, leaf_list){
  l_sel <- length(sel)
  if(l_sel== 0){ return(NA)}
  if(type_borne == "HMM"){
    return( min(which(a >= 1 - alpha) -1, l_sel))
  }
  
  if(type_borne == "HMM_small"){
    
    return( max(c(which(a < alpha),1)) - 1)
  }
  if(type_borne == "Simes"){
    # pS <- posthocBySimes(pval, sel, alpha, Rcpp = FALSE, verbose = FALSE) 
    #
    k_val <- sapply(1:l_sel, function(k){ sum(pval[sel] > alpha * k / m) + k - 1})
    return( min(k_val))
    # return(length(sel)- pS)
  }
  if(type_borne=="DKW_tree"){
    return( V.star(sel, C, ZL, leaf_list))
    
    
  }
  if(type_borne=="DKW"){
    # return( V.star(sel, C, ZL, leaf_list))
    C <- sqrt(0.5 * log(1 / alpha))
    pval_sel <- pval[sel]
    pi <- pval_sel[order(pval_sel)]
    DKW_fun <- function(i) {
      (C / (2 * (1 - pi[i])) +
         (C^2 / (4 * (1 - pi[i])^2) +
            (l_sel - i) / (1 - pi[i]))^(1 / 2))^2
    }
    return(min(c(sapply(1:l_sel, DKW_fun), l_sel)))
    
  }
}
