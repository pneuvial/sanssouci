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
