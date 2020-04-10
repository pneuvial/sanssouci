#' Title
#'
#' @param m
#' @param alpha
#' @param beta
#' @param A
#' @param Pi
#' @param f0x
#' @param f1x
#'
#' @return
#' @export
#'
#' @examples
monte_carl_borne <- function(m, alpha, beta, A, Pi, f0x, f1x, B) {
  borne_list <- future_map(1:B, ~ sim_x_kn(m, alpha, beta, A, Pi, f0x, f1x))
  return(map_dbl(borne_list, ~ sum(. == 0)))
}
