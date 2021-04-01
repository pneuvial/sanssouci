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
boots_param_unknown_f0 <- function(A_est,  Pi_est, x_from, prob1, h, Sel_from, al,
                                   seuil, min_size,  min_jump = NULL, n,  max_pi0, m0_init, sd0_init, df_init, 
                                   norm_init, type_init, approx){
  Sel_from <- Sel_from %>% rename(Sel_from = Sel)
  m <- length(x_from)
  Data_temp <- sim_hmm_from_weightkde( A_est,  Pi_est,  x_from, prob1, h, n )
  theta <- Data_temp$theta 
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
  rm(Em)
  gc()
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
  
  Sel <- Selection_tibble(x, fw_bc_EM_star, seuil, A_est_star, f0x_est_star, f1x_est_star,
                          Pi_est_sar, min_size, min_jump) %>% 
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
boots_param_known_f0 <- function (A_est, Pi_est, x_from, prob1, h, Sel_from, al, seuil, 
                                  min_size, n, max_pi0, m0_init, sd0_init, df_init, norm_init, 
                                  type_init, approx, maxit=100) 
{
  Sel_from <- Sel_from %>% rename(Sel_from = Sel)
  m <- length(x_from)
  Data_temp <- sim_hmm_from_weightkde(A_est, Pi_est, x_from, 
                                      prob1, h, n)
  theta <- Data_temp$theta
  x <- Data_temp$x
  rm(Data_temp)
  gc()
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
    rm(f1x_first)
    gc()
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
           x, eps = 1e-04, maxit = maxit, h = h, f0_known = TRUE, 
           approx)
  A_est_star <- Em$A
  Pi_est_sar <- Em$Pi
  f1x_est_star <- Em$f1x
  f0x_est_star <- f0x
  fw_bc_EM_star <- Em$fw_bc_EM
  rm(Em)
  gc()
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
    rm(f1x_first)
    gc()
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
  Sel <- Selection_tibble(x, fw_bc_EM_star, seuil, A_est_star, 
                          f0x_est_star, f1x_est_star, Pi_est_sar, min_size,
                          pval = pval) %>% 
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


