
#' Title
#'
#' @param x 
#' @param p0 
#'
#' @return
#' @export
#'
#' @examples
pval_m <-function(x, p0){
  ord  <- order(x)
  CDF1 <- cumsum(p0[ord]) / sum(p0[ord])
  CDF  <- cumsum(p0[ord]) / sum(p0[ord])
  CDF[CDF > 0.5] <- 1 - CDF[CDF > 0.5]
  pvalm <- 2 * CDF
  return(pvalm[order(ord)])
}


#' Title
#'
#' @param x 
#' @param fw_bc 
#' @param seuil 
#' @param A_est 
#' @param f0x_est 
#' @param f1x_est 
#' @param Pi_est 
#' @param min_size 
#'
#' @return
#' @export
#'
#' @examples
Selection_tibble <- function(x, fw_bc, seuil, A_est, f0x_est, f1x_est, Pi_est, min_size){
  
  m <- length(x)
  
  LIS <- enframe(fw_bc$gamma[, 1]) %>%
    arrange(value) %>%
    rowid_to_column() %>%
    mutate(k = cumsum(value) / rowid) %>%
    filter(k <= seuil) %>%
    arrange(name)
  
  sel_sc <- LIS$name
  viterbi <- viterbi_log(m, log(A_est), log(f0x_est), 
                         log(f1x_est), log(Pi_est))
  
  
  sel_viter_min_size <- long_reg(viterbi, min_size)
  sel_viter_est <- which(viterbi == 1)
  pvalm <- pval_m(x, fw_bc$gamma[,1])
  pvalm_tresh <- which(pvalm < seuil)
  
  Sel <- tibble( Sel = list(
    pvalm_tresh,
    sel_sc, 
    sel_viter_est,
    sel_viter_min_size),
    Nom = c("pval_tresh",
            "lfdr_tresh", 
            "sel_viter_est",
            "sel_viter_min_size"
    )) %>%
    mutate(Size = map_dbl(Sel,~length(.)))
}