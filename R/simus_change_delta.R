
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
#' @param min_jump 
#' @param pval 
#' @param all 
#'
#' @return
#' @export
#'
#' @examples
Selection_delta <- function(x, fw_bc, seuil, A_est,
                            f0x_est, f1x_est, Pi_est, min_size,
                            min_jump = NULL,
                            pval = NULL, all = FALSE){

  m <- length(x)

  LIS <- enframe(fw_bc$gamma[, 1]) %>%
    arrange(value) %>%
    rowid_to_column() %>%
    mutate(k = cumsum(value) / rowid) %>%
    filter(k <= seuil) %>%
    arrange(name)

  sel_sc <- LIS$name
  

  sel_viter_est <- which(viterbi == 1)
  if(is.null(pval)){
    pval <- pval_m(x, fw_bc$gamma[,1])
  }
  pvalm_tresh <- which(pval < seuil)
  pvalm_sel <- pvalm_tresh[pvalm_tresh < 200]


  Sel <- tibble( Sel = list(
    pvalm_sel,
    pvalm_tresh,
    sel_sc),
    Nom = c(
      "pvalm_sel",
      "pval_tresh",
      "lfdr_tresh")) %>%
    mutate(Size = map_dbl(Sel,~length(.)))
}


#' Title
#'
#' @param A_est
#' @param Pi_est
#' @param x_from
#' @param prob1
#' @param h
#' @param Sel_from
#' @param al
#' @param seuil
#' @param min_size
#' @param n
#' @param max_pi0
#' @param m0_init
#' @param sd0_init
#' @param df_init
#' @param norm_init
#' @param type_init
#' @param approx
#' @param maxit
#' @param f0_known
#' @param delta
#'
#' @return
#' @export
#'
#' @examples
#' m <- 2000
#' theta <- sim_markov(m, Pi = c(0.8,0.2), A = matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T))
#' x <- rep(0, m)
#' x[theta == 0] <- rnorm(sum(theta ==0))
#' x[theta == 1] <- rnorm(sum(theta ==1), 2, 1)
#' Est <- Estimation(x, m0_init = 0, sd0_init = 1, norm_init = TRUE, plot= TRUE)
#' Sel <- Selection_delta(x, Est$Em$fw_bc_EM, 0.05,
#'  Est$Em$A,  Est$Em$f0x,  Est$Em$f1x,
#'   Est$Em$Pi, 9)
#' boots_delta(Est$Em$A, Est$Em$Pi,
#' x_from = x,
#' prob1 = Est$Em$fw_bc_EM$gamma[,2], h,
#' Sel, al =0.1,
#' seuil=0.05,
#' min_size =9,b_act*s_dbnr,
#' n = n,
#' max_pi0 = 0.99999,
#' m0_init =  0, sd0_init = 1,
#'  norm_init = TRUE,
#' type_init = "given",
#' approx = TRUE,
#' delta = 0.9)
boots_delta <- function (A_est, Pi_est, x_from, prob1, h, Sel_from, al, seuil,
                         min_size, n, max_pi0, m0_init, sd0_init, df_init, norm_init,
                         type_init, approx, maxit=100,f0_known = TRUE,delta)
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
    if(!f0_known){
      f0x_first <- sapply(d0$x, function(xi) {
        sum(K((x_from - xi)/h) * (1 - prob1))/sum(h * (1 -
                                                         prob1))
      })
      f0x <- approx(d1$x, f0x_first, x)$y
      rm(f0x_first)
    }
    gc()
  }
  else {
    f1x <- sapply(x, function(xi) {
      sum(K((x_from - xi)/h) * prob1)/sum(h * prob1)
    })
    if(!f0_known){
      f0x <- sapply(x, function(xi) {
        sum(K((x_from - xi)/h) * prob1)/sum(h * (1-prob1))
      })
    }
  }
  fw_bc_or_star <- for_back(m, A_est, f0x, f1x, Pi_est)
  Pis_or_star <- lapply(2:m, function(i) {
    get_A(m, alpha = fw_bc_or_star$alpha, beta = fw_bc_or_star$beta,
          A_est, f0x, f1x, i = i)
  })

  Est <- Estimation( x, h =h,
                     m0_init, sd0_init, df_init, norm_init, max_pi0= max_pi0,
                     f0_known = f0_known, f0x_est = NULL, pval = NULL,
                     plot = FALSE, size_plot= min(10000, length(x)),
                     approx = approx)

  if(f0_known){
    pval_don <- pval
  }else{
    pval_don = NULL
  }
  A_est <- Est$Em$A
  Pi_est <- Est$Em$Pi
  Sel <- Selection_delta(x, Est$Em$fw_bc_EM, seuil, A_est,
                         f0x_est =Est$Em$f0x ,
                         f1x_est= Est$Em$f1x, Pi_est, min_size,
                         pval = pval_don,
                         all = all ) %>%
    rename(Size_boot = Size)
  Pis_est_star <- lapply(2:m, function(i) {
    get_A(m, alpha = Est$fw_bc_EM$alpha, beta = Est$fw_bc_EM$beta,
          A_est, Est$Em$f0x, Est$Em$f1x, i = i)
  })
  if (approx) {
    d1_from <- density(x_from)
    f1x_first <- sapply(d1_from$x, function(xi) {
      sum(K((x - xi)/h) * Est$fw_bc_EM$gamma[, 2])/sum(h * Est$fw_bc_EM$gamma[, 2])
    })
    f1x_from <- approx(d1_from$x, f1x_first, x_from)$y
    rm(f1x_first)
    gc()
  }
  else {
    f1x_from <- sapply(x_from, function(xi) {
      sum(K((x - xi)/h) * Est$fw_bc_EM$gamma[, 2])/sum(h * Est$fw_bc_EM$gamma[, 2])
    })
    f0x_from <- sapply(x_from, function(xi) {
      sum(K((x - xi)/h) * Est$fw_bc_EM$gamma[, 1])/sum(h * Est$fw_bc_EM$gamma[, 1])
    })
  }
  fw_bc_from <- for_back(m, A_est, f0x_from, f1x_from,
                         Pi_est)
  Pis_from <- lapply(2:m, function(i) {
    get_A(m, alpha = fw_bc_from$alpha, beta = fw_bc_from$beta,
          A_est, f0x_from, f1x_from, i = i)
  })
  prob <- c(al/2, 1-al/2,
            al*(1-delta), 1 - al*(1-delta),
            al, 1-al)
  Sel_boot <- Sel %>% full_join(Sel_from, by = "Nom") %>%
    mutate(Borne_from = map(Sel_from,
                            ~get_probs(sel = ., li0 = fw_bc_from$gamma[, 1],
                                       Pis = Pis_from, f0x = f0x_from, f1x = f1x_from, prob)),
           Borne_oracle_boot = map(Sel, ~get_probs(sel = .,
                                                   li0 = fw_bc_or_star$gamma[, 1], Pis = Pis_or_star,
                                                   f0x = f0x, f1x = f1x, prob)),
           Borne_est_boot = map(Sel,   ~get_probs(sel = ., li0 = fw_bc_EM_star$gamma[, 1],
                                                  Pis = Pis_est_star, f0x = f0x_est_star,
                                                  f1x = f1x_est_star, prob)),
           V_HMM_oracle_boot = map_dbl(Borne_oracle_boot, ~.[6]),
           V_HMM_small_oracle_boot = map_dbl(Borne_oracle_boot, ~.[5]),
           V_HMM_oracle_boot_aldemi = map_dbl(Borne_oracle_boot, ~.[2]),
           V_HMM_small_oracle_boot_aldemi = map_dbl(Borne_oracle_boot, ~.[1]),
           V_HMM_oracle_boot_aldelta = map_dbl(Borne_oracle_boot, ~.[4]),
           V_HMM_small_oracle_boot_aldelta = map_dbl(Borne_oracle_boot, ~.[3]),
           V_HMM_est_boot_aldemi = map_dbl(Borne_est_boot, ~.[2]),
           V_HMM_small_est_boot_aldemi = map_dbl(Borne_est_boot, ~.[1]),
           V_HMM_est_boot_aldelta = map_dbl(Borne_est_boot, ~.[4]),
           V_HMM_small_est_boot_aldelta = map_dbl(Borne_est_boot, ~.[3]),
           V_HMM_est_boot = map_dbl(Borne_est_boot, ~.[6]),
           V_HMM_small_est_boot = map_dbl(Borne_est_boot, ~.[5]),
           V_HMM_est_boot_aldemi_samesel = map_dbl(Borne_from, ~.[2]),
           V_HMM_small_est_boot_aldemi_samesel = map_dbl(Borne_from, ~.[1]),
           V_HMM_est_boot_aldelta_samesel = map_dbl(Borne_from, ~.[4]),
           V_HMM_small_est_boot_aldelta_samesel = map_dbl(Borne_from, ~.[3]),
           V_HMM_est_boot_samesel = map_dbl(Borne_from, ~.[6]),
           V_HMM_small_est_boot_samesel = map_dbl(Borne_from, ~.[5]),
           Espe = map_dbl(Sel, ~sum(fw_bc_or_star$gamma[.,1])),
           Real_boot = map_dbl(Sel, ~sum(theta[.] == 0))) %>%
    select(-Borne_oracle_boot, -Borne_est_boot, -Borne_from)
  return(Sel_boot)
}


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
#' @param type_init
#' @param num_seed
#' @param f0_known
#' @param approx
#' @param all
#' @param size_b0
#' @param pct_b1
#' @param include_H0
#' @param delta
#'
#' @return
#' @export
#'
#' @examples
#' simu_delta( m = c(100),
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
#' norm_init = TRUE, df= 2, num_seed= 1234, type_init="given", f0_known=TRUE,
#' approx = TRUE, delta = 0.9)
simu_delta <- function(m, A, Pi, n, rho, SNR, prob, type_sim = "HMM", al, s_dbnr,
                       b_act, d, seuil,
                       h =0.3,  n_boot, min_size, norm, m0, sd0, df,
                       m0_init, sd0_init, df_init, norm_init, max_pi0= 0.99999,
                       type_init, num_seed, f0_known, approx, all =FALSE,
                       size_b0= 300, pct_b1 =1/3, include_H0 = FALSE, delta) {
  set.seed(num_seed)
  ## Simuation des données
  if(type_sim =="HMM"){
    theta <- sim_markov(m, Pi, A)
  }
  if(type_sim =="block"){
    theta <- rep(rep(0:1, c(size_b0, size_b0 * pct_b1)),
                 m / (size_b0 * (1 + pct_b1)))
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

  if(f0_known){
    pval_don <- pval
  }else{
    pval_don = NULL
  }
  A_est <- Est$Em$A
  Pi_est <- Est$Em$Pi
  H0 <- sum(theta == 0)
  H1 <- sum(theta == 1)
  Sel <- Selection_delta(x, Est$Em$fw_bc_EM, seuil, A_est,
                         f0x_est =Est$Em$f0x ,
                         f1x_est= Est$Em$f1x, Pi_est, min_size,
                         pval = pval_don,
                         all = all )

  Pis_est <- lapply(2:m, function(i){
    get_A( m,alpha = Est$Em$fw_bc_EM$alpha, beta = Est$Em$fw_bc_EM$beta,
           A_est, f0x =Est$Em$f0x ,
           f1x= Est$Em$f1x, i = i)
  })


  ### Pour Estimer bootstrap :
  boots <- enframe( x = 1:n_boot, name = NULL, value = "id_boot") %>%
    mutate(HMM_boot = map(id_boot, ~boots_delta(A_est, Pi_est,
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
                                                approx = approx,
                                                delta = delta))) %>%
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
      Borne_est = map(Sel,~get_probs(sel = ., li0 = Est$Em$fw_bc_EM$gamma[,1],
                                     Pis = Pis_est, f0x =Est$Em$f0x ,
                                     f1x= Est$Em$f1x, probs = c(al/2, 1-al/2,
                                                                al, 1-al,
                                                                al*(1-delta), 1 - al*(1-delta)))),
      V_simes = map_dbl(Sel,~borne ( type_borne = "Simes", sel = ., m = m,
                                     pval = pval, alpha = al)
      ),
      V_HMM_small_or =  map_dbl(IC_or,~.[1]),
      V_HMM_or =  map_dbl(IC_or,~.[2]),
      V_HMM_est_aldemi =  map_dbl(Borne_est,~.[2]),
      V_HMM_small_est_aldemi =  map_dbl(Borne_est,~.[1]),
      V_HMM_est_aldelta =  map_dbl(Borne_est,~.[6]),
      V_HMM_small_est_aldelta =  map_dbl(Borne_est,~.[5]),
      V_HMM_small_est =  map_dbl(Borne_est,~.[3]),
      V_HMM_est =  map_dbl(Borne_est,~.[4]),
      FDR_est = map_dbl(Sel, ~sum( Est$Em$fw_bc_EM$gamma[.,1])),
      FDR_or  = map_dbl(Sel, ~sum(fw_bc_or$gamma[.,1])),
      V_real = map_dbl(Sel , ~sum(theta[.] == 0))
    ) %>%
    select(-Sel,  -Borne_est) %>%
    left_join(boots, by = c("Nom"))

  return(Final)
}





# tic()
#
# Res_simu_delta <- Param_easy %>%
#   mutate(borne = future_pmap(
#     list(m = m,  Pi = Pi, A = A, n = n, rho = rho,
#          SNR = SNR, prob = prob, type_sim = type_sim, al = al,
#          n_boot = n_boot,
#          s_dbnr = s_dbnr, b_act = b_act, d =d, seuil =seuil,
#          min_size = min_size,  h =0.3,
#          norm = norm, m0 = m0, sd0 =sd0, df = df,
#          sd0_init = sd0_init, m0_init= m0_init, norm_init = norm_init,
#          num_seed = simu*3+5, type_init="given", approx = TRUE,
#          f0_known= f0_known, delta = delta),
#     simu_delta
#   )) %>%
#   unnest(borne)
# toc()
#
# save(Res_simu_delta, file = "res_simu_delta.RData")

