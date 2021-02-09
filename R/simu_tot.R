
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

simu_tot <- function(m, A, Pi, n, rho, SNR, prob, type_sim, al, s_dbnr,
                     b_act, d, seuil,
                     h =0.3,  n_boot, min_size, norm, m0, sd0, df,
                     m0_init, sd0_init, df_init, norm_init, max_pi0= 0.99999, 
                     type_init, num_seed, f0_known, approx) {
  set.seed(num_seed)
  theta <- sim_markov(m,Pi, A)
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
  
  if(norm_init){
    f0x_est <- dnorm(x, m0_init, sd0_init)
  }else{
    f0x_est <- dt(x, df_init)
  }
  pi0_hat <- max(min(sum(pval > 0.8) / (m * 0.2), max_pi0), 0.6)
  f_hatx <- x %>%
    map_dbl( ~f_hatK(x, ., h = h,K) )
  f1x_est <-  f1x_hat(f0x_est, f_hatx, pi0_hat)
  f1x_est[f1x_est <= 0] <- min(f0x_est)
  mini <- max(0.6, ((1 + max_pi0) * pi0_hat -max_pi0) / pi0_hat)
  a <-runif(1, mini, max_pi0)
  b <- 1 - a
  c <- pi0_hat * b / (1 - pi0_hat)
  d <- 1 - c
  Em <- Em(m, A = matrix(c(a, b, c, d), byrow = TRUE, ncol=2),
                  Pi= c(pi0_hat, 1 -pi0_hat),  f0x = f0x_est, f1x = f1x_est,
                  x, eps = 0.0001,
                  maxit =1000, h = h, f0_known = f0_known, approx = approx)
  # if(Em$A[1,1] > Em$A[2,2]){ 
    
    A_est <- Em$A
    Pi_est <-Em$Pi
    f0x_est <- Em$f0x
    f1x_est <- Em$f1x
    fw_bc_EM <- Em$fw_bc_EM
  # }else{ 
    A_est <- Em$A[2:1,2:1]
    Pi_est <-Em$Pi[2:1]
    f0x_est <- Em$f1x
    f1x_est <- Em$f0x
    fw_bc_EM <- for_back(m, A_est,f0x_est, f1x_est, Pi_est)
  # }
  
  Pis_est <- lapply(2:m, function(i){
    get_A( m,alpha = fw_bc_EM$alpha, beta = fw_bc_EM$beta, A_est, f0x_est, 
           f1x_est, i = i)
  })
  if(f0_known){
    boots_param <- boots_param_known_f0
  }else{
    boots_param <- boots_param_unknown_f0
  }
 Sel <- Selection_tibble(x, fw_bc_EM, seuil, A_est, f0x_est, f1x_est, Pi_est, min_size)
   ### Pour Estimer bootstrap : 
  boots <- enframe( x = 1:n_boot, name = NULL, value = "id_boot") %>% 
    mutate(HMM_boot = map(id_boot, ~boots_param(A_est, Pi_est, 
                                                x_from = x,
                                                prob1 = fw_bc_EM$gamma[,2], h, 
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
      IC_or_aldemi = map(Sel,~get_IC(sel = ., li0 = fw_bc_or$gamma[,1],
                                      Pis = Pis_or, f0x =f0x ,
                                      f1x= f1x, alpha = al/2)),
      IC_or = map(Sel,~get_IC(sel = ., li0 = fw_bc_or$gamma[,1],
                              Pis = Pis_or, f0x =f0x ,
                              f1x= f1x, alpha = al)),
      IC_est_aldemi = map(Sel,~get_IC(sel = ., li0 = fw_bc_EM$gamma[,1],
                                      Pis = Pis_est, f0x =f0x_est ,
                                      f1x= f1x_est, alpha = al/2)),
      IC_est = map(Sel,~get_IC(sel = ., li0 = fw_bc_EM$gamma[,1],
                               Pis = Pis_est, f0x =f0x_est ,
                               f1x= f1x_est, alpha = al)),
      V_simes = map_dbl(Sel,~borne ( type_borne = "Simes", sel = ., m = m,
                                     pval = pval, alpha = al)
      ),
      V_HMM_or_aldemi =  map_dbl(IC_or_aldemi,~.[2]),
      V_HMM_small_or_aldemi =  map_dbl(IC_or_aldemi,~.[1]),
      V_HMM_small_or =  map_dbl(IC_or,~.[1]),
      V_HMM_or =  map_dbl(IC_or,~.[2]),
      V_HMM_est_aldemi =  map_dbl(IC_est_aldemi,~.[2]),
      V_HMM_small_est_aldemi =  map_dbl(IC_est_aldemi,~.[1]),
      V_HMM_small_est =  map_dbl(IC_est,~.[1]),
      V_HMM_est =  map_dbl(IC_est,~.[2]),
      FDR_est = map_dbl(Sel, ~sum( fw_bc_EM$gamma[.,1])),
      FDR_or  = map_dbl(Sel, ~sum(fw_bc_or$gamma[.,1])),
      V_real = map_dbl(Sel , ~sum(theta[.] == 0))
    ) %>%
    select(-Sel,  -IC_est, -IC_est_aldemi) %>% 
    left_join(boots, by = c("Nom"))
  
  return(Final)
}