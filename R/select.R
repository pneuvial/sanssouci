#' Title
#'
#' @param type_sel 
#' @param EM 
#' @param stuff 
#' @param min_size 
#'
#' @return
#' @export
#'
#' @examples
selection <- function(type_sel, EM, stuff, seuil_lis, seuil_p, min_size = 0 ){
  if(type_sel =="lis"){
  LIS <- enframe(EM$fw_bc_EM$gamma[, 1]) %>%
    arrange(value) %>%
    rowid_to_column() %>%
    mutate(k = cumsum(value) / rowid) %>%
    filter(k <= seuil_lis) %>%
    arrange(name)
return(LIS$name)
    }
  if(type_sel == "viter"){
    viterbi <- viterbi_log(stuff$m, log(EM$A), log(stuff$f0x),
                           log(EM$f1x), log(EM$Pi)) 
    # long_reg(viterbi, min_size)
    return(which(viterbi == 1) )
    }  
 if(type_sel=="pval"){
   which(stuff$pval < seuil_p)  
 }
}