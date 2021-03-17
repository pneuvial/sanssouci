
Selection_add_knowledge <- function(x, fw_bc, seuil, A_est,
         f0x_est, f1x_est, Pi_est, min_size, 
         min_jump = NULL,
         pval = NULL, all = FALSE, H0){
  
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
  if(!is.null(min_jump)){
    diff <- sel_viter_min_size[-1] - sel_viter_min_size[-length(sel_viter_min_size)]
    tomerge <- which(diff < min_jump & round(diff)!=1)
    if(length(tomerge)>=1){
      sel_viter_min_size_supp <- unlist(lapply(tomerge,function(x){
        (sel_viter_min_size[x] +1): (sel_viter_min_size[x+1]-1 )
      } ))
      sel_viter_min_size <- sort(c(sel_viter_min_size,sel_viter_min_size_supp))
    }
    
  }
  sel_viter_est <- which(viterbi == 1)
  if(is.null(pval)){
    pval <- pval_m(x, fw_bc$gamma[,1])
  }
  pvalm_tresh <- which(pval < seuil)
  sel_viter_est_H0 <- which(viterbi == 1 & H0)
  pvalm_tresh_H0 <- which(pval < seuil & H0)
  H1 <- wich(!H0)
  if(all){
    Sel <- tibble( Sel = list(
      1:m,
      pvalm_tresh,
      sel_sc, 
      sel_viter_est,
      sel_viter_min_size),
      Nom = c("all",
              "pval_tresh",
              "lfdr_tresh", 
              "sel_viter_est",
              "sel_viter_min_size"
      )) %>%
      mutate(Size = map_dbl(Sel,~length(.)))
  }else{
    
    Sel <- tibble( Sel = list(
      pvalm_tresh,
      sel_sc, 
      sel_viter_est,
      sel_viter_min_size,
      sel_viter_est_H0,
      pvalm_tresh_H0, H1),
      Nom = c(
        "pval_tresh",
        "lfdr_tresh", 
        "sel_viter_est",
        "sel_viter_min_size",
        "sel_viter_est_H0",
        "pvalm_tresh_H0","H1"
      )) %>%
      mutate(Size = map_dbl(Sel,~length(.)))
  }
}
## Bootstrap


#' Title
#'
#' @param A_est the  transition matrix
#' @param Pi_est initial probabilities (if null the stationnary distributio  is used)
#' @param x_from the observed value
#' @param prob1 the probability of being in state one ( obtain using forward backward algorithm for instance)
#' @param h the size of the window in the kernel algorithm
#' @param Sel_from the selection tibble for theobserved value (output of selection_tibble )
#' @param al the given risk
#' @param seuil threshold for the pvalues
#' @param min_size the minimum size of concomitant selected position in S ( by default ony used for viterbi_min_size)
#' @param n the number of individuals
#' @param max_pi0 a given maximum values of pi_0
#' @param m0_init if the "guess" distribution of x under the nulll state is the normal this is the guess value of the mean
#' @param sd0_init if the "guess" distribution of x under the nulll state is the normal this is the guess value of the standard
#' @param df_init if the distribution given to the EM is student, this is the degreee of freedom
#' @param norm_init logical, is the distribution under the null given to the EM algo is gaussian or not ? (if not consider as student)
#'
#' @return
#' @export
#'
#' @examples
knowledge_boots_param_unknown_f0 <- function(A_est,  Pi_est, x_from, prob1, h, Sel_from, al,
                                   seuil, min_size,  min_jump = NULL, n,  max_pi0, m0_init, sd0_init, df_init, 
                                   norm_init, type_init, approx){
  Sel_from <- Sel_from %>% rename(Sel_from = Sel)
  m <- length(x_from)
  Data_temp <- sim_hmm_from_weightkde( A_est,  Pi_est,  x_from, prob1, h, n )
  theta <- Data_temp$theta 
  H0 <- (theta==0)
  x <- Data_temp$x  
  if(norm_init){
    pval <- 2 * (1 - pnorm(abs(x), m0_init, sd0_init))
  }else{
    pval <- 2 * (1 - pt(abs(x), df_init)) 
  }
  
  
  if(type_init == "locfdr"){
    w <- locfdr(x)
    m0_init <- w$fp0[3,1]
    sd0_init <-  w$fp0[3,2]
    norm_init <- TRUE 
  }
  
  if (approx) {
    d0 <- density(x)
    d1 <- density(x)
    f0x_first <- sapply(d0$x, function(xi) {
      sum(K((x_from - xi)/h) * (1 - prob1))/sum(h * (1 - 
                                                       prob1))
    })
    f1x_first <- sapply(d1$x, function(xi) {
      sum(K((x_from - xi)/h) * prob1)/sum(h * prob1)
    })
    f1x <- approx(d1$x, f1x_first, x)$y
    f0x <- approx(d0$x, f0x_first, x)$y
  }else{
    f0x  <-sapply(x, function(xi){
      sum(K((x_from -xi)/h) * (1-prob1)) / sum(h *(1- prob1))
    })
    f1x  <-sapply(x, function(xi){
      sum(K((x_from -xi)/h) * prob1) / sum(h * prob1)
    }) 
  }
  
  
  
  ## oracle boot 
  fw_bc_or_star <- for_back(m, A_est, f0x, f1x, Pi_est)
  Pis_or_star <- lapply(2:m, function(i){
    get_A( m,alpha = fw_bc_or_star$alpha, beta = fw_bc_or_star$beta,
           A_est, f0x, f1x, i = i)
  } )
  
  ## est boot
  if(norm_init){
    f0x_est_star <- dnorm(x, m0_init, sd0_init)
  }else{
    f0x_est_star <- dt(x, df_init)
  }
  pi0_hat <- max(min(sum(pval > 0.8) / (m * 0.2), max_pi0 ), 0.6)
  if(approx){
    d <- density(x,bw = h)
    f_hatx <- approx(d$x,d$y,x)$y
  }else {
    f_hatx <- x %>%
      map_dbl( ~f_hatK(x, ., h = h,K) )
  }
  f1x_est_star <-  f1x_hat(f0x_est_star, f_hatx, pi0_hat)
  f1x_est_star[f1x_est_star <= 0] <- min(f0x_est_star)
  mini <- max(0.6, ((1 + max_pi0) * pi0_hat -max_pi0) / pi0_hat)
  a <-runif(1, mini, max_pi0)
  b <- 1 - a
  c <- pi0_hat * b / (1 - pi0_hat)
  d <- 1 - c
  
  if(max(c(a,b,c,d))> 1 ){ stop("Pb de A init (boot)")}
  
  Em <- Em(m, A = matrix(c(a, b, c, d), byrow = TRUE, ncol=2),
           Pi= c(pi0_hat, 1 -pi0_hat),  f0x = f0x_est_star, f1x = f1x_est_star,
           x, eps = 0.0001,
           maxit =1000, h = h, f0_known = FALSE, approx)
  # if(Em$A[1,1] > Em$A[2,2]){ 
  A_est_star <- Em$A
  Pi_est_sar <-Em$Pi
  f1x_est_star <- Em$f1x
  f0x_est_star <- Em$f0x
  fw_bc_EM_star <- Em$fw_bc_EM
  # }else{ 
  #   A_est_star <- Em$A[2:1,2:1]
  #   Pi_est_sar <-Em$Pi[2:1]
  #   f0x_est_star <- Em$f1x
  #   f1x_est_star <- Em$f0x
  #   fw_bc_EM_star <- for_back(m, A_est_star,f0x_est_star, f1x_est_star, Pi_est_sar)
  # }
  
  Pis_est_star <- lapply(2:m, function(i){
    get_A( m,alpha = fw_bc_EM_star$alpha, beta = fw_bc_EM_star$beta, 
           A_est_star, f0x_est_star, f1x_est_star, i = i)
  })
  
  gamma_EM_star1 <- fw_bc_EM_star$gamma[,2]
  gamma_EM_star0 <- fw_bc_EM_star$gamma[,1]
  
  if (approx) {
    d0_from <- density(x_from)
    d1_from <- density(x_from)
    f0x_first <- sapply(d0_from$x, function(xi) {
      sum(K((x - xi)/h) * gamma_EM_star0)/sum(h * gamma_EM_star0)
    })
    f1x_first <- sapply(d1_from$x, function(xi) {
      sum(K((x - xi)/h) * gamma_EM_star1)/sum(h * gamma_EM_star1)
    })
    f1x_from <- approx(d1_from$x, f1x_first, x_from)$y
    f0x_from <- approx(d0_from$x, f0x_first, x_from)$y
  }else{
    f1x_from  <-sapply(x_from, function(xi){
      sum(K((x -xi)/h) * gamma_EM_star1) / sum(h * gamma_EM_star1)
    })
    f0x_from  <-sapply(x_from, function(xi){
      sum(K((x -xi)/h) * gamma_EM_star0) / sum(h * gamma_EM_star0)
    })  
  }
  
  fw_bc_from <-  for_back(m, A_est_star, f0x_from, f1x_from, Pi_est_sar)
  Pis_from <- lapply(2:m, function(i){
    get_A( m,alpha = fw_bc_from$alpha, beta = fw_bc_from$beta, 
           A_est_star, f0x_from, f1x_from, i = i)
  })
  ## selection SC 
  
  Sel <- Selection_add_knowledge(x, fw_bc_EM_star, seuil, A_est_star, f0x_est_star, f1x_est_star,
                          Pi_est_sar, min_size, min_jump, H0) %>% 
    rename(Size_boot = Size)
  
  # Sel_boot <- Sel %>% left_join(Sel_from,  by = "Nom") %>% 
  #   mutate(
  #     quantile_from = map(Sel_from,~get_quantiles(
  #       sel = ., li0 =fw_bc_from$gamma[,1],
  #       Pis = Pis_from, f0x = f0x_from, f1x = f1x_from)),
  #     quantile_oracle_boot = map(Sel,~get_quantiles(
  #       sel = ., li0 =fw_bc_or_star$gamma[,1],
  #       Pis = Pis_or_star, f0x = f0x, f1x = f1x)), 
  #     quantile_est_boot = map(Sel,~get_quantiles(
  #       sel = ., li0 =fw_bc_EM_star$gamma[,1],
  #       Pis = Pis_est_star, f0x = f0x_est_star, f1x = f1x_est_star
  #     )),
  #     V_HMM_oracle_boot_aldemi =  map2_dbl(Sel, quantile_oracle_boot,~borne(
  #       type_borne = "HMM", sel = .x, a = .y, alpha = al /2)
  #     ),
  #     V_HMM_small_oracle_boot_aldemi =  map2_dbl(Sel, quantile_oracle_boot,~borne(
  #       type_borne = "HMM_small", sel = .x, a = .y, alpha = al / 2)
  #     ),
  #     V_HMM_est_boot_aldemi =  map2_dbl(Sel, quantile_est_boot,~borne(
  #       type_borne = "HMM", sel = .x, a = .y, alpha = al /2)
  #     ),
  #     V_HMM_small_est_boot_aldemi =  map2_dbl(Sel, quantile_est_boot,~borne(
  #       type_borne = "HMM_small", sel = .x, a = .y, alpha = al /2)
  #     ),
  #     V_HMM_est_boot_aldemi_samesel =  map2_dbl(Sel_from, quantile_from,~borne(
  #       type_borne = "HMM", sel = .x, a = .y, alpha = al /2)
  #     ),
  #     V_HMM_small_est_boot_aldemi_samesel =  map2_dbl(Sel_from, quantile_from,~borne(
  #       type_borne = "HMM_small", sel = .x, a = .y, alpha = al /2)
  #     ),
  #     Espe = map_dbl(Sel, ~sum(fw_bc_or_star$gamma[.,1])),
  #     Real_boot = map_dbl(Sel , ~sum(theta[.] == 0))
  #   ) %>% select(-quantile_oracle_boot, -quantile_est_boot, -quantile_from)
  
  
  Sel_boot <- Sel %>% full_join(Sel_from, by = "Nom") %>%
    mutate(IC_from = map(Sel_from, 
                         ~get_IC(sel = ., li0 = fw_bc_from$gamma[, 1], 
                                 Pis = Pis_from, f0x = f0x_from, f1x = f1x_from, al/2)), 
           IC_oracle_boot = map(Sel, ~get_IC(sel = ., 
                                             li0 = fw_bc_or_star$gamma[, 1], Pis = Pis_or_star, 
                                             f0x = f0x, f1x = f1x, al/2)), 
           IC_est_boot = map(Sel,   ~get_IC(sel = ., li0 = fw_bc_EM_star$gamma[, 1], 
                                            Pis = Pis_est_star, f0x = f0x_est_star, f1x = f1x_est_star, al/2)), 
           V_HMM_oracle_boot_aldemi = map_dbl(IC_oracle_boot, ~.[2]), 
           V_HMM_small_oracle_boot_aldemi = map_dbl(IC_oracle_boot, ~.[1]),
           V_HMM_est_boot_aldemi = map_dbl(IC_est_boot, ~.[2]), 
           V_HMM_small_est_boot_aldemi = map_dbl(IC_est_boot, ~.[1]), 
           V_HMM_est_boot_aldemi_samesel = map_dbl(IC_from, ~.[2]),
           V_HMM_small_est_boot_aldemi_samesel = map_dbl(IC_from, ~.[1]),
           Espe = map_dbl(Sel, ~sum(fw_bc_or_star$gamma[.,1])), 
           Real_boot = map_dbl(Sel, ~sum(theta[.] == 0))) %>% 
    select(-IC_oracle_boot, -IC_est_boot, -IC_from)
  return(Sel_boot)
  
}


## Bootstrap


#' Title
#'
#' @param A_est the  transition matrix
#' @param Pi_est initial probabilities (if null the stationnary distributio  is used)
#' @param x_from the observed value
#' @param prob1 the probability of being in state one ( obtain using forward backward algorithm for instance)
#' @param h the size of the window in the kernel algorithm
#' @param Sel_from the selection tibble for theobserved value (output of selection_tibble )
#' @param al the given risk
#' @param seuil threshold for the pvalues
#' @param min_size the minimum size of concomitant selected position in S ( by default ony used for viterbi_min_size)
#' @param n the number of individuals
#' @param max_pi0 a given maximum values of pi_0
#' @param m0_init if the "guess" distribution of x under the nulll state is the normal this is the guess value of the mean
#' @param sd0_init if the "guess" distribution of x under the nulll state is the normal this is the guess value of the standard
#' @param df_init if the distribution given to the EM is student, this is the degreee of freedom
#' @param norm_init logical, is the distribution under the null given to the EM algo is gaussian or not ? (if not consider as student)
#'
#' @return
#' @export
#'
#' @examples
knowledge_boots_param_known_f0 <- function (A_est, Pi_est, x_from, prob1, h, Sel_from, al, seuil, 
                                  min_size, n, max_pi0, m0_init, sd0_init, df_init, norm_init, 
                                  type_init, approx, H0) 
{
  Sel_from <- Sel_from %>% rename(Sel_from = Sel)
  m <- length(x_from)
  Data_temp <- sim_hmm_from_weightkde(A_est, Pi_est, x_from, 
                                      prob1, h, n)
  theta <- Data_temp$theta
  H0 <- (theta==0)
  de
  x <- Data_temp$x
  if (norm_init) {
    pval <- 2 * (1 - pnorm(abs(x), m0_init, sd0_init))
  }
  else {
    pval <- 2 * (1 - pt(abs(x), df_init))
  }
  if (type_init == "locfdr") {
    w <- locfdr(x)
    m0_init <- w$fp0[3, 1]
    sd0_init <- w$fp0[3, 2]
    norm_init <- TRUE
  }
  if (norm_init) {
    f0x <- dnorm(x, m0_init, sd0_init)
    f0x_from <- dnorm(x_from, m0_init, sd0_init)
  }
  else {
    f0x <- dt(x, df_init)
    f0x_from <- dt(x_from, df_init)
  }
  
  if (approx) {
    d1 <- density(x)
    
    f1x_first <- sapply(d1$x, function(xi) {
      sum(K((x_from - xi)/h) * prob1)/sum(h * prob1)
    })
    f1x <- approx(d1$x, f1x_first, x)$y
  }
  else {
    f1x <- sapply(x, function(xi) {
      sum(K((x_from - xi)/h) * prob1)/sum(h * prob1)
    })
  }
  fw_bc_or_star <- for_back(m, A_est, f0x, f1x, Pi_est)
  Pis_or_star <- lapply(2:m, function(i) {
    get_A(m, alpha = fw_bc_or_star$alpha, beta = fw_bc_or_star$beta, 
          A_est, f0x, f1x, i = i)
  })
  f0x_est_star <- f0x
  pi0_hat <- max(min(sum(pval > 0.8)/(m * 0.2), max_pi0), 0.6)
  if (approx) {
    d <- density(x, bw = h)
    f_hatx <- approx(d$x, d$y, x)$y
  }
  else {
    f_hatx <- x %>% map_dbl(~f_hatK(x, ., h = h, K))
  }
  f1x_est_star <- f1x_hat(f0x_est_star, f_hatx, pi0_hat)
  f1x_est_star[f1x_est_star <= 0] <- min(f0x_est_star)
  mini <- max(0.6, ((1 + max_pi0) * pi0_hat - max_pi0)/pi0_hat)
  a <- runif(1, mini, max_pi0)
  b <- 1 - a
  c <- pi0_hat * b/(1 - pi0_hat)
  d <- 1 - c
  if (max(c(a, b, c, d)) > 1) {
    stop("Pb de A init (boot)")
  }
  Em <- Em(m, A = matrix(c(a, b, c, d), byrow = TRUE, ncol = 2), 
           Pi = c(pi0_hat, 1 - pi0_hat), f0x = f0x_est_star, 
           f1x = f1x_est_star, 
           x, eps = 1e-04, maxit = 1000, h = h, f0_known = TRUE, 
           approx)
  A_est_star <- Em$A
  Pi_est_sar <- Em$Pi
  f1x_est_star <- Em$f1x
  f0x_est_star <- f0x
  fw_bc_EM_star <- Em$fw_bc_EM
  Pis_est_star <- lapply(2:m, function(i) {
    get_A(m, alpha = fw_bc_EM_star$alpha, beta = fw_bc_EM_star$beta, 
          A_est_star, f0x_est_star, f1x_est_star, i = i)
  })
  gamma_EM_star1 <- fw_bc_EM_star$gamma[, 2]
  gamma_EM_star0 <- fw_bc_EM_star$gamma[, 1]
  if (approx) {
    d1_from <- density(x_from)
    f1x_first <- sapply(d1_from$x, function(xi) {
      sum(K((x - xi)/h) * gamma_EM_star1)/sum(h * gamma_EM_star1)
    })
    f1x_from <- approx(d1_from$x, f1x_first, x_from)$y
  }
  else {
    f1x_from <- sapply(x_from, function(xi) {
      sum(K((x - xi)/h) * gamma_EM_star1)/sum(h * gamma_EM_star1)
    })
    f0x_from <- sapply(x_from, function(xi) {
      sum(K((x - xi)/h) * gamma_EM_star0)/sum(h * gamma_EM_star0)
    })
  }
  fw_bc_from <- for_back(m, A_est_star, f0x_from, f1x_from, 
                         Pi_est_sar)
  Pis_from <- lapply(2:m, function(i) {
    get_A(m, alpha = fw_bc_from$alpha, beta = fw_bc_from$beta, 
          A_est_star, f0x_from, f1x_from, i = i)
  })
  Sel <- Selection_add_knowledge(x, fw_bc_EM_star, seuil, A_est_star, 
                          f0x_est_star, f1x_est_star, Pi_est_sar, min_size,H0) %>% 
    rename(Size_boot = Size)
  Sel_boot <- Sel %>% full_join(Sel_from, by = "Nom") %>%
    mutate(IC_from = map(Sel_from, 
                         ~get_IC(sel = ., li0 = fw_bc_from$gamma[, 1], 
                                 Pis = Pis_from, f0x = f0x_from, f1x = f1x_from, al/2)), 
           IC_oracle_boot = map(Sel, ~get_IC(sel = ., 
                                             li0 = fw_bc_or_star$gamma[, 1], Pis = Pis_or_star, 
                                             f0x = f0x, f1x = f1x, al/2)), 
           IC_est_boot = map(Sel,   ~get_IC(sel = ., li0 = fw_bc_EM_star$gamma[, 1], 
                                            Pis = Pis_est_star, f0x = f0x_est_star, f1x = f1x_est_star, al/2)), 
           V_HMM_oracle_boot_aldemi = map_dbl(IC_oracle_boot, ~.[2]), 
           V_HMM_small_oracle_boot_aldemi = map_dbl(IC_oracle_boot, ~.[1]),
           V_HMM_est_boot_aldemi = map_dbl(IC_est_boot, ~.[2]), 
           V_HMM_small_est_boot_aldemi = map_dbl(IC_est_boot, ~.[1]), 
           V_HMM_est_boot_aldemi_samesel = map_dbl(IC_from, ~.[2]),
           V_HMM_small_est_boot_aldemi_samesel = map_dbl(IC_from, ~.[1]),
           Espe = map_dbl(Sel, ~sum(fw_bc_or_star$gamma[.,1])), 
           Real_boot = map_dbl(Sel, ~sum(theta[.] == 0))) %>% 
    select(-IC_oracle_boot, -IC_est_boot, -IC_from)
  return(Sel_boot)
}


# quantile_from = map(Sel_from,~get_quantiles(
#   sel = ., li0 =fw_bc_from$gamma[,1],
#   Pis = Pis_from, f0x = f0x_from, f1x = f1x_from)),
# quantile_oracle_boot = map(Sel,~get_quantiles(
#   sel = ., li0 =fw_bc_or_star$gamma[,1],
#   Pis = Pis_or_star, f0x = f0x, f1x = f1x)), 
# quantile_est_boot = map(Sel,~get_quantiles(
#   sel = ., li0 =fw_bc_EM_star$gamma[,1],
#   Pis = Pis_est_star, f0x = f0x_est_star, f1x = f1x_est_star
# )),

#' Title
#'
#' @param m 
#' @param A 
#' @param Pi 
#' @param n 
#' @param rho 
#' @param SNR 
#' @param prob 
#' @param type_sim 
#' @param al 
#' @param s_dbnr 
#' @param b_act 
#' @param d 
#' @param seuil 
#' @param h 
#' @param n_boot 
#' @param min_size 
#' @param norm 
#' @param m0 
#' @param sd0 
#' @param df 
#' @param m0_init 
#' @param sd0_init 
#' @param df_init 
#' @param norm_init 
#' @param max_pi0 
#'
#' @return
#' @export
#'
#' @examples
#' simu_tot( m = c(100),
#' A = matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T),
#' Pi = c(0.95, 0.05),
#' n = c(500),
#' rho = c(0),
#' SNR = 2,
#' prob = c(0.5),
#' type_sim = c("HMM"),
#' n_boot = 20,
#' al = 0.2, s_dbnr = 10, b_act = 2, d = 1, seuil= 0.05,
#' min_size = 2, norm = TRUE, sd0 = 0.5, m0= 0,sd0_init = 0.5, m0_init= 0,
#'  norm_init = TRUE, df= 2, num_seed= 1234, type_init="given", f0_known=FALSE, 
#'  approx = TRUE)

simu_tot_H0 <- function(m, A, Pi, n, rho, SNR, prob, type_sim = "HMM", al, s_dbnr,
                     b_act, d, seuil,
                     h =0.3,  n_boot, min_size, norm, m0, sd0, df,
                     m0_init, sd0_init, df_init, norm_init, max_pi0= 0.99999, 
                     type_init, num_seed, f0_known, approx, all =FALSE,
                     size_b0= 300, pct_b1 =1/3, include_H0 = FALSE) {
  set.seed(num_seed)
  if(type_sim =="HMM"){
    theta <- sim_markov(m, Pi, A)
  }
  if(type_sim =="block"){
    theta <- rep(rep(0:1, c(size_b0, size_b0*pct_b1)),m/(size_b0*(1+pct_b1))) 
    nb0 <-  m / (1+pct_b1)
    nb1 <- m - nb0
    A <- matrix(c((nb0-8)/(nb0-1), (8)/(nb0-1),(7)/(nb1-1),(nb1-7)/(nb1-1)), ncol =2, byrow = TRUE)
  }
  x <- rep(0, m)
  if(norm){
    x[theta == 0] <- rnorm(sum(theta ==0), m0, sd0)
    x[theta == 1] <- rnorm(sum(theta ==1), SNR*sd0 + m0, sd0)
  }else{
    x[theta == 0] <- rt(sum(theta ==0),df)
    x[theta == 1] <- rt(sum(theta ==1),df) + SNR
  }
  
  if(type_init == "locfdr"){
    w <- locfdr(x)
    m0_init <- w$fp0[3,1]
    sd0_init <-  w$fp0[3,2]
    norm_init <- TRUE 
  }
  
  if(norm_init){
    pval <- 2 * (1 - pnorm(abs(x), m0_init, sd0_init))
  }else{
    pval <- 2 * (1 - pt(abs(x), df_init)) 
  }
  
  
  ## Pour HMM oracle
  if(norm){
    f0x <- dnorm(x, m0, sd0)
    f1x <- dnorm(x, SNR*sd0 + m0, sd0)
  }else{
    f0x <- dt(x, df)
    f1x <- dt(x-SNR, df)
  }
  
  fw_bc_or <- for_back(m, A, f0x, f1x, Pi)
  Pis_or <- lapply(2:m, function(i){
    get_A( m,alpha = fw_bc_or$alpha, beta = fw_bc_or$beta, A, f0x, f1x, i = i)
  } )
  
  
  ## Pour HMM est
  
  Est <- Estimation( x, h =h,
                     m0_init, sd0_init, df_init, norm_init, max_pi0= max_pi0, 
                     f0_known = f0_known, f0x_est = NULL, pval = NULL, 
                     plot = FALSE, size_plot= min(10000, length(x)), 
                     approx = approx)
  
  
  A_est <- Est$Em$A
  Pi_est <- Est$Em$Pi
 H0 <- (theta==0)
  Sel <- Selection_add_knowledge(x, Est$Em$fw_bc_EM, seuil, A_est, 
                          f0x_est =Est$Em$f0x , 
                          f1x_est= Est$Em$f1x, Pi_est, min_size, all = all ,H0) 
  
  Pis_est <- lapply(2:m, function(i){
    get_A( m,alpha = Est$Em$fw_bc_EM$alpha, beta = Est$Em$fw_bc_EM$beta,
           A_est, f0x =Est$Em$f0x ,
           f1x= Est$Em$f1x, i = i)
  })
  
  # 
  # Estimation (x,  h =0.3,
  #                        m0_init, sd0_init, df_init, norm_init, max_pi0, 
  #                        f0_known = TRUE, f0x_est = NULL, pval = NULL, 
  #                        plot = FALSE, size_plot= min(10000, length(x)), 
  #                        approx = TRUE)
  # if(norm_init){
  #   f0x_est <- dnorm(x, m0_init, sd0_init)
  # }else{
  #   f0x_est <- dt(x, df_init)
  # }
  # pi0_hat <- max(min(sum(pval > 0.8) / (m * 0.2), max_pi0), 0.6)
  # # f_hatx <- x %>%
  # #   map_dbl( ~f_hatK(x, ., h = h,K) ) ## 
  # 
  # if(approx){
  #   d <- density(x,bw = h)
  #   f_hatx <- approx(d$x,d$y,x)$y
  # }else {
  #   f_hatx <- x %>%
  #     map_dbl( ~f_hatK(x, ., h = h,K) )
  # }
  # 
  # f1x_est <-  f1x_hat(f0x_est, f_hatx, pi0_hat)
  # f1x_est[f1x_est <= 0] <- min(f0x_est)
  # mini <- max(0.6, ((1 + max_pi0) * pi0_hat -max_pi0) / pi0_hat)
  # a <-runif(1, mini, max_pi0)
  # b <- 1 - a
  # c <- pi0_hat * b / (1 - pi0_hat)
  # d <- 1 - c
  # Em <- Em(m, A = matrix(c(a, b, c, d), byrow = TRUE, ncol=2),
  #                 Pi= c(pi0_hat, 1 -pi0_hat),  f0x = f0x_est, f1x = f1x_est,
  #                 x, eps = 0.0001,
  #                 maxit =1000, h = h, f0_known = f0_known, approx = approx)
  # # if(Em$A[1,1] > Em$A[2,2]){ 
  #   
  #   A_est <- Em$A
  #   Pi_est <-Em$Pi
  #   f0x_est <- Em$f0x
  #   f1x_est <- Em$f1x
  #   fw_bc_EM <- Em$fw_bc_EM
  # # }else{ 
  #   A_est <- Em$A[2:1,2:1]
  #   Pi_est <-Em$Pi[2:1]
  #   f0x_est <- Em$f1x
  #   f1x_est <- Em$f0x
  #   fw_bc_EM <- for_back(m, A_est,f0x_est, f1x_est, Pi_est)
  # # }
  # 
  # Pis_est <- lapply(2:m, function(i){
  #   get_A( m,alpha = fw_bc_EM$alpha, beta = fw_bc_EM$beta, A_est, f0x_est, 
  #          f1x_est, i = i)
  # })
  if(f0_known){
    boots_param <- knowledge_boots_param_known_f0
  }else{
    boots_param <- boots_param_unknown_f0
  }
  ### Pour Estimer bootstrap : 
  boots <- enframe( x = 1:n_boot, name = NULL, value = "id_boot") %>% 
    mutate(HMM_boot = map(id_boot, ~boots_param(A_est, Pi_est, 
                                                x_from = x,
                                                prob1 = Est$Em$fw_bc_EM$gamma[,2], h, 
                                                Sel, al, 
                                                seuil, 
                                                min_size,b_act*s_dbnr, 
                                                n = n, 
                                                max_pi0 = max_pi0,
                                                m0_init =  m0_init, sd0_init = sd0_init,
                                                df_init  = df_init, norm_init = norm_init, 
                                                type_init = type_init,
                                                approx = approx))) %>% 
    unnest(HMM_boot) %>% 
    select(- Sel) %>% 
    nest(Est_HMM_boot = c(id_boot, Real_boot,
                          V_HMM_est_boot_aldemi,
                          V_HMM_small_est_boot_aldemi,
                          V_HMM_oracle_boot_aldemi,
                          V_HMM_small_oracle_boot_aldemi,
                          V_HMM_est_boot_aldemi_samesel,
                          V_HMM_small_est_boot_aldemi_samesel,
                          Espe, 
                          Size, 
                          Size_boot, 
                          Sel_from)) 
  
  Det_a <- det(A_est)
  
  Final <- Sel   %>%
    mutate(
      sd0_init_est = sd0_init, 
      m0_init_est = m0_init,
      Det_A_est = Det_a,
      IC_or = map(Sel,~get_IC(sel = ., li0 = fw_bc_or$gamma[,1],
                              Pis = Pis_or, f0x =f0x ,
                              f1x= f1x, alpha = al)),
      IC_est_aldemi = map(Sel,~get_IC(sel = ., li0 = Est$Em$fw_bc_EM$gamma[,1],
                                      Pis = Pis_est, f0x =Est$Em$f0x ,
                                      f1x= Est$Em$f1x, alpha = al/2)),
      IC_est = map(Sel,~get_IC(sel = ., li0 = Est$Em$fw_bc_EM$gamma[,1],
                               Pis = Pis_est, f0x =Est$Em$f0x ,
                               f1x= Est$Em$f1x, alpha = al)),
      V_simes = map_dbl(Sel,~borne ( type_borne = "Simes", sel = ., m = m,
                                     pval = pval, alpha = al)
      ),
      V_HMM_small_or =  map_dbl(IC_or,~.[1]),
      V_HMM_or =  map_dbl(IC_or,~.[2]),
      V_HMM_est_aldemi =  map_dbl(IC_est_aldemi,~.[2]),
      V_HMM_small_est_aldemi =  map_dbl(IC_est_aldemi,~.[1]),
      V_HMM_small_est =  map_dbl(IC_est,~.[1]),
      V_HMM_est =  map_dbl(IC_est,~.[2]),
      FDR_est = map_dbl(Sel, ~sum( Est$Em$fw_bc_EM$gamma[.,1])),
      FDR_or  = map_dbl(Sel, ~sum(fw_bc_or$gamma[.,1])),
      V_real = map_dbl(Sel , ~sum(theta[.] == 0))
    ) %>%
    select(-Sel,  -IC_est, -IC_est_aldemi) %>% 
    left_join(boots, by = c("Nom"))
  
  return(Final)
}
