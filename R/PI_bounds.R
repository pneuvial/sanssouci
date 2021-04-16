#' Title
#'
#' @param Sel 
#' @param Est 
#' @param al 
#' @param delta 
#'
#' @return
#' @export
#'
#' @examples
PI_bounds <- function(Sel, Est, al, delta,Pis_est , pval){
  probs <- c(al/2, 1-al/2,
    al, 1-al,
    al*(1-delta), 1 - al*(1-delta), 
    al*(delta), 1 - al*(delta))
  m <- length( pval)
  Final <- Sel   %>%
    mutate(
      Det_A_est = det(Est$Em$A),
      Borne_est = map(Sel,~get_probs(sel = ., li0 = Est$Em$fw_bc_EM$gamma[,1],
                                     Pis = Pis_est, f0x =Est$Em$f0x ,
                                     f1x= Est$Em$f1x, probs = probs)),
      V_simes = map_dbl(Sel,~borne ( type_borne = "Simes", sel = ., m = m,
                                     pval = pval, alpha = al)
      ),
      V_HMM_est_aldemi =  map_dbl(Borne_est,~.[2]),
      V_HMM_small_est_aldemi =  map_dbl(Borne_est,~.[1]),
      V_HMM_est_aldelta =  map_dbl(Borne_est,~.[7]),
      V_HMM_small_est_aldelta =  map_dbl(Borne_est,~.[8]),
      V_HMM_est_al1_moins_delta =  map_dbl(Borne_est,~.[6]),
      V_HMM_small_est_al1_moins_delta =  map_dbl(Borne_est,~.[5]),
      V_HMM_small_est =  map_dbl(Borne_est,~.[3]),
      V_HMM_est =  map_dbl(Borne_est,~.[4]),
      FDR_est = map_dbl(Sel, ~sum( Est$Em$fw_bc_EM$gamma[.,1]))
    ) %>%
    select(-Sel,  -Borne_est) 
    return(Final)
}