#' Title
#'
#' @param A_est matrix. The estimated transition matrix
#' @param f0 numeric vector of size 2 (the number of states). The mean and sd of the null hypothesis (supposed to be gaussian)
#' @param Pi_est numeric vector of size 2 (the number of states). The estimated initial probabilities
#' @param x_from numeric vector. The observed variables
#' @param prob1 numeric vector. The estimated probability for each observation to be alternatives. 
#' @param h numeric. Size of the window for the kernel estimation of densities
#' @param Sel_from tibble. The different type of selection.
#' @param al numeric, the threshold of 
#' @param seuil threshold of selected p-values
#' @param min_size the minimum size 
#' @param num_b 
#' @param n 
#' @param f0x_from 
#' @param type_sel 
#'
#' @return
#' @export
#'
#' @examples
boots <- function(A_est, f0, Pi_est, x_from, prob1, h, Sel_from, al,
                        seuil, min_size,  f0x_from,type_sel,num_b=2, n=2){
  Sel_from <- Sel_from %>% rename(Sel_from = Sel)
  m <- length(x_from)
  Data_temp <- sim_hmm_from_weightkde( A_est, f0, Pi_est,  x_from, prob1, h, n )
  theta <- Data_temp$theta 
  x <- Data_temp$x
  pval <- 2 * (1 - pnorm(abs(x)))
  f0x <- dnorm(x, f0[1], f0[2])
  f1x  <-sapply(x, function(xi){
    sum(K((x_from -xi)/h) * prob1) / sum(h * prob1)
  })
  ## oracle boot 
  fw_bc_or_star <- for_back(m, A_est, f0x, f1x, Pi_est)
  Pis_or_star <- lapply(2:m, function(i){
    get_A( m,alpha = fw_bc_or_star$alpha, beta = fw_bc_or_star$beta,
           A_est, f0x, f1x, i = i)
  } )
  
  ## est boot
  pi0_hat <- min(sum(pval > 0.8) / (m * 0.2), 0.97)
  f_hatx <- x %>%
    map_dbl( ~f_hatK(x, ., h = h,K) )
  f1x_est_star <-  f1x_hat(f0x, f_hatx, pi0_hat)
  f1x_est_star[f1x_est_star <= 0] <- min(f0x)
  
  a <- runif(1, 0.7, 0.98)
  b <- 1 - a
  c <- runif(1, 0.1, 0.3)
  d <- 1 - c
  
  Em <- Em_tot(m, A = matrix(c(a, b, c, d), byrow = TRUE, ncol=2),
               Pi= c(pi0_hat, 1 -pi0_hat),  f0x = f0x, f1x = f1x_est_star,
               x, eps = 0.0001,
               maxit =1000, h = h)
  
  A_est_star <- Em$A
  Pi_est_sar <-Em$Pi[,1]
  f1x_est_star <- Em$f1x
  fw_bc_EM_star <- Em$fw_bc_EM
  Pis_est_star <- lapply(2:m, function(i){
    get_A( m,alpha = fw_bc_EM_star$alpha, beta = fw_bc_EM_star$beta, 
           A_est_star, f0x, f1x_est_star, i = i)
  })
  
  gamma_EM_star <- fw_bc_EM_star$gamma[,2]
  
  f1x_from  <-sapply(x_from, function(xi){
    sum(K((x -xi)/h) * gamma_EM_star) / sum(h * gamma_EM_star)
  })
  fw_bc_from <-  for_back(m, A_est_star, f0x_from, f1x_from, Pi_est_sar)
  Pis_from <- lapply(2:m, function(i){
    get_A( m,alpha = fw_bc_from$alpha, beta = fw_bc_from$beta, 
           A_est_star, f0x_from, f1x_from, i = i)
  })
  ## selection SC 
  LIS <- enframe(fw_bc_EM_star$gamma[, 1]) %>%
    arrange(value) %>%
    rowid_to_column() %>%
    mutate(k = cumsum(value) / rowid) %>%
    filter(k <= seuil) %>%
    arrange(name)
  
  sel_sc <- LIS$name
  
  seuil_pval <- sum(fw_bc_EM_star$gamma[pval < seuil, 1]) / sum(pval < seuil)
  LIS_idpval <- enframe(fw_bc_EM_star$gamma[, 1]) %>%
    arrange(value) %>%
    rowid_to_column() %>%
    mutate(k = cumsum(value) / rowid) %>%
    filter(k <= seuil_pval) %>%
    arrange(name)
  
  sel_idpval <- LIS_idpval$name
  
  
  
  viterbi <- viterbi_log(m, log(A_est_star), log(f0x), 
                         log(f1x_est_star), log(Pi_est))
  
  seuil_viter <- sum(fw_bc_EM_star$gamma[viterbi == 1,1]) / sum(viterbi ==1)
  
  LIS_viter <- enframe(fw_bc_EM_star$gamma[, 1]) %>%
    arrange(value) %>%
    rowid_to_column() %>%
    mutate(k = cumsum(value) / rowid) %>%
    filter(k <= seuil_viter) %>%
    arrange(name)
  
  sel_lis_viter <- LIS_viter$name
  sel_viter_min_size <- long_reg(viterbi, min_size)
  pval01 <- pval
  pval01[pval01 < seuil] <- 1
  pval01[pval01 > seuil] <- 0
  
  
  sel_viter_est <- which(viterbi == 1)
  H0 <- which(theta == 0)
  H1 <- which(theta == 1)
  pval_tresh <- which(pval < seuil)
  pval_lfdr <- sort(order(pval)[1:length(sel_sc)])
  
  Sel <- tibble( Sel = list(
    H1,
    pval_tresh,
    pval_lfdr, 
    sel_sc, 
    sel_idpval,
    sel_lis_viter,
    sel_viter_est,
    sel_viter_min_size,
    intersect(H0, pval_tresh),
    intersect(H0, sel_viter_est),
    c(1:num_b)),
    Type_sel = c("H1",
            "pval_tresh",
            "pval_ord",
            "lfdr_tresh", 
            "lfdr_tresh_pval",
            "lfdr_tresh_viter",
            "sel_viter_est",
            "sel_viter_min_size",
            "H0_and_pval_tresh_small",
            "H0_and_sel_viter_est",
            "block")) %>%
    mutate( Size_boot = map_dbl(Sel,~length(.))) %>% 
    filter(Type_sel %in% type_sel)
  
  
  Sel_boot <- Sel %>% left_join(Sel_from,  by = "Type_sel") %>% 
    mutate(
      quantile_from = map(Sel_from,~get_quantiles(
        sel = ., li0 =fw_bc_from$gamma[,1],
        Pis = Pis_from, f0x = f0x_from, f1x = f1x_from)),
      quantile_oracle_boot = map(Sel,~get_quantiles(
        sel = ., li0 =fw_bc_or_star$gamma[,1],
        Pis = Pis_or_star, f0x = f0x, f1x = f1x)), 
      quantile_est_boot = map(Sel,~get_quantiles(
        sel = ., li0 =fw_bc_EM_star$gamma[,1],
        Pis = Pis_est_star, f0x = f0x, f1x = f1x_est_star
      )),
      V_HMM_oracle_boot_aldemi =  map2_dbl(Sel, quantile_oracle_boot,~borne(
        type_borne = "HMM", sel = .x, a = .y, alpha = al /2)
      ),
      V_HMM_small_oracle_boot_aldemi =  map2_dbl(Sel, quantile_oracle_boot,~borne(
        type_borne = "HMM_small", sel = .x, a = .y, alpha = al / 2)
      ),
      V_HMM_est_boot_aldemi =  map2_dbl(Sel, quantile_est_boot,~borne(
        type_borne = "HMM", sel = .x, a = .y, alpha = al /2)
      ),
      V_HMM_small_est_boot_aldemi =  map2_dbl(Sel, quantile_est_boot,~borne(
        type_borne = "HMM_small", sel = .x, a = .y, alpha = al /2)
      ),
      V_HMM_est_boot_aldemi_samesel =  map2_dbl(Sel_from, quantile_from,~borne(
        type_borne = "HMM", sel = .x, a = .y, alpha = al /2)
      ),
      V_HMM_small_est_boot_aldemi_samesel =  map2_dbl(Sel_from, quantile_from,~borne(
        type_borne = "HMM_small", sel = .x, a = .y, alpha = al /2)
      ),
      Espe = map_dbl(Sel, ~sum(fw_bc_or_star$gamma[.,1])),
      Real_boot = map_dbl(Sel , ~sum(theta[.] == 0))
    ) %>% select(-quantile_oracle_boot, -quantile_est_boot, -quantile_from)
  return(Sel_boot)
  
}