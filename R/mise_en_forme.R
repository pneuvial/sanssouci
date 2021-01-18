#' Title
#'
#' @param Result_hmm_boot 
#'
#' @return
#' @export
#'
#' @examples
res_all_stat <- function(Result_hmm_boot){
  
  alpha <- Result_hmm_boot$al %>% unique()
  Res_all <-   Result_hmm_boot %>% 
    filter(Nom != "pval_tresh_min_size" & Size > 0) %>% 
    mutate(Low_boot = map_dbl(Est_HMM_boot, function(x){
      mutate(x, eps =   Real_boot - Espe) %>% 
        pull(eps) %>% 
        quantile(probs=alpha, type = 3)
    }),
    High_boot = map_dbl(Est_HMM_boot, function(x){
      mutate(x, eps =   Real_boot - Espe) %>% 
        pull(eps) %>% 
        quantile(probs=(1-alpha), type = 3)
    }),
    q_aldemi_high = map_dbl(Est_HMM_boot, function(x){
      mutate(x, eps =   (V_HMM_oracle_boot_aldemi - V_HMM_est_boot_aldemi) ) %>% 
        pull(eps) %>% 
        quantile(probs=(1 - alpha / 2), na.rm =TRUE, type = 3)
    }),
    q_aldemi_low = map_dbl(Est_HMM_boot, function(x){
      mutate(x, eps =   V_HMM_small_oracle_boot_aldemi - V_HMM_small_est_boot_aldemi) %>% 
        pull(eps) %>% 
        quantile(probs= (alpha / 2), na.rm =TRUE, type = 3)
    }),
    q_samesel = map2_dbl(Est_HMM_boot, V_HMM_est_aldemi, function(x,y){
      mutate(x, eps =  y -  V_HMM_est_boot_aldemi_samesel ) %>% 
        pull(eps) %>% 
        quantile(probs= (1- alpha / 2), na.rm =TRUE, type = 3)
    }),
    q_samesel_small = map2_dbl(Est_HMM_boot, V_HMM_small_est_aldemi, function(x,y){
      mutate(x, eps =  y-  V_HMM_small_est_boot_aldemi_samesel ) %>% 
        pull(eps) %>% 
        quantile(probs= (alpha / 2), na.rm =TRUE, type = 3)
    }),
    V_HMM_boot_samesel = V_HMM_est_aldemi + q_samesel,
    V_HMM_small_boot_samesel = V_HMM_small_est_aldemi + q_samesel_small,
    V_HMM_boot_compReal = FDR_est + High_boot,
    V_HMM_small_boot_compReal = FDR_est + Low_boot, 
    V_HMM_boot_qdemi = V_HMM_est_aldemi + q_aldemi_high,
    V_HMM_small_boot_qdemi = V_HMM_small_est_aldemi + q_aldemi_low, 
    V_HMM_boot_naif = map_dbl(Est_HMM_boot, ~ pull(. , Real_boot) %>% 
                                quantile(probs = (1 - alpha), type = 3)),
    V_HMM_small_boot_naif = map_dbl(Est_HMM_boot, ~ pull(. , Real_boot) %>% 
                                      quantile(probs = (alpha), type = 3)), 
    FDR_boot_naif = map_dbl(Est_HMM_boot, ~ pull(. , Real_boot) %>% 
                              mean()),
    Nom = case_when( Nom == "pval_tresh" ~  "S(X) == bgroup('{',p < 0.05,'}')", 
                     Nom == "sel_viter_est" ~ "S(X) == bgroup('{',viterbi == 1,'}')",
                     Nom == "sel_viter_min_size" ~ "S(X) == bgroup('{','viterbi = 1 & size > 4','}')",
                     Nom == "H0_and_pval_tresh_small" ~  "S(X) ==H0~intersect()~ bgroup('{',p < 0.05 ,'}')", 
                     Nom == "H0_and_sel_viter_est" ~ "S(X) ==H0~intersect()~bgroup('{',viterbi == 1 ,'}')",
                     Nom == "H1" ~  "S(X) == H1", 
                     Nom == "pval_ord" ~  "S(X) == bgroup('{',p < th,'}')",
                     Nom == "lfdr_tresh" ~  "S(X) == SC(0.05)",
                     Nom == "lfdr_tresh_pval" ~  "S(X) == SC(FDR[p])",
                     Nom == "lfdr_tresh_viter"~  "S(X) == SC(FDR[v])",
                     Nom == "block" ~ "S(X) ==bgroup('{','1, ..., 200' ,'}')"), 
    Nom = factor(Nom, levels = c("S(X) == bgroup('{',p < 0.05,'}')","S(X) == SC(FDR[p])",
                                 "S(X) == SC(0.05)", "S(X) == bgroup('{',p < th,'}')",
                                 "S(X) == bgroup('{',viterbi == 1,'}')", 
                                 "S(X) == SC(FDR[v])",
                                 "S(X) == bgroup('{','viterbi = 1 & size > 4','}')", 
                                 "S(X) ==bgroup('{','1, ..., 200' ,'}')",
                                 "S(X) == H1",
                                 "S(X) ==H0~intersect()~ bgroup('{',p < 0.05 ,'}')", 
                                 "S(X) ==H0~intersect()~bgroup('{',viterbi == 1 ,'}')"))
    ) 
  return(Res_all) 
}

#' Title
#'
#' @param Result_hmm_boot 
#'
#' @return
#' @export
#'
#' @examples
res_all_FDR <- function(Result_hmm_boot){
  alpha <- Result_hmm_boot$al %>% unique()
  Res_all <-   Result_hmm_boot %>% 
    filter(Nom != "pval_tresh_min_size" & Size > 0) %>% 
    mutate(Low_boot = map_dbl(Est_HMM_boot, function(x){
      vec <- mutate(x, eps =  ( Real_boot - Espe)/Size_boot) %>% 
        pull(eps) %>% sort()
      vec[floor(length(vec) * alpha)]
      # quantile(probs=alpha, type = 3, na.rm =TRUE)
    }),
    High_boot = map_dbl(Est_HMM_boot, function(x){
      mutate(x, eps =   (Real_boot - Espe)/Size_boot) %>% 
        pull(eps) %>% 
        quantile(probs=(1-alpha), type = 3, na.rm =TRUE)
    }),
    q_aldemi_high = map_dbl(Est_HMM_boot, function(x){
      mutate(x, eps =   (V_HMM_oracle_boot_aldemi - V_HMM_est_boot_aldemi)/Size_boot ) %>% 
        pull(eps) %>% 
        quantile(probs=(1 - alpha / 2), na.rm =TRUE, type = 3)
    }),
    q_aldemi_low = map_dbl(Est_HMM_boot, function(x){
      vec <- x %>% filter(Size_boot >0) %>%
        mutate( eps =   (V_HMM_small_oracle_boot_aldemi - V_HMM_small_est_boot_aldemi)/Size_boot) %>% 
        pull(eps) %>% sort()
      
      vec[max(floor(length(vec) * alpha/2),1)]
      # quantile(probs= (alpha / 2), na.rm =TRUE, type = 3)
    }),
    q_samesel = map2_dbl(Est_HMM_boot, V_HMM_est_aldemi, function(x,y){
      mutate(x, eps =  y -  V_HMM_est_boot_aldemi_samesel ) %>% 
        pull(eps) %>% 
        quantile(probs= (1- alpha / 2), na.rm =TRUE, type = 3)
    }),
    q_samesel_small = map2_dbl(Est_HMM_boot, V_HMM_small_est_aldemi, function(x,y){
      vec<-x %>% filter(Size_boot >0) %>%
        mutate(  eps =  y-  V_HMM_small_est_boot_aldemi_samesel ) %>% 
        pull(eps) %>% sort() 
      
      vec[max(floor(length(vec) * alpha/2),1)]
      # quantile(probs= (alpha / 2), na.rm =TRUE, type = 3)
    }),
    V_HMM_boot_samesel = V_HMM_est_aldemi + q_samesel,
    V_HMM_boot_samesel= V_HMM_boot_samesel /Size,
    V_HMM_small_boot_samesel = V_HMM_small_est_aldemi + q_samesel_small,
    V_HMM_small_boot_samesel= V_HMM_small_boot_samesel/Size,
    V_HMM_boot_compReal = FDR_est/Size + High_boot,
    V_HMM_small_boot_compReal = FDR_est/Size + Low_boot, 
    V_HMM_boot_qdemi = V_HMM_est_aldemi/Size + q_aldemi_high,
    V_HMM_small_boot_qdemi = V_HMM_small_est_aldemi/Size + q_aldemi_low, 
    V_HMM_boot_naif = map_dbl(Est_HMM_boot, ~pull(. , Real_boot) %>% 
                                quantile(probs = (1 - alpha), type = 3, na.rm =TRUE)),
    V_HMM_boot_naif= V_HMM_boot_naif/Size,
    V_HMM_small_boot_naif = map_dbl(Est_HMM_boot, function(x){
      vec <-  pull(x , Real_boot) %>% sort()
      vec[floor(length(vec) * alpha)]
    }),
    # quantile(probs = (alpha), type = 3, na.rm =TRUE)),
    V_HMM_small_boot_naif = V_HMM_small_boot_naif/Size,
    FDR_boot_naif = map_dbl(Est_HMM_boot, ~ pull(. , Real_boot) %>% 
                              mean()),
    V_HMM_est =V_HMM_est / Size, 
    V_HMM_oracle = V_HMM_oracle / Size, 
    V_HMM_semi_oracle = V_HMM_semi_oracle/ Size, 
    V_HMM_small_est =V_HMM_small_est / Size, 
    V_HMM_small_oracle = V_HMM_small_oracle / Size, 
    V_HMM_small_semi_oracle = V_HMM_small_semi_oracle/ Size,
    V_DKW_tree = V_DKW_tree /Size, 
    V_simes = V_simes / Size, 
    V_real = V_real / Size, 
    FDR_boot_naif = FDR_boot_naif /Size, 
    FDR_or = FDR_or/ Size , FDR_est = FDR_est/Size,
    Nom = case_when( Nom == "pval_tresh" ~  "S(X) == bgroup('{',p < 0.05,'}')", 
                     Nom == "sel_viter_est" ~ "S(X) == bgroup('{',viterbi == 1,'}')",
                     Nom == "sel_viter_or" ~ "S(X) == bgroup('{',viterbi[or] == 1,'}')",
                     Nom == "sel_viter_min_size" ~ "S(X) == bgroup('{','viterbi = 1 & size > 4','}')",
                     Nom == "H0_and_pval_tresh_small" ~  "S(X) ==H0~intersect()~ bgroup('{',p < 0.05 ,'}')", 
                     Nom == "H0_and_sel_viter_est" ~ "S(X) ==H0~intersect()~bgroup('{',viterbi == 1 ,'}')",
                     Nom == "H1" ~  "S(X) == H1", 
                     Nom == "pval_ord" ~  "S(X) == bgroup('{',p < th,'}')",
                     Nom == "lfdr_tresh" ~  "S(X) == SC(0.05)",
                     Nom == "lfdr_tresh_pval" ~  "S(X) == SC(FDR[p])",
                     Nom == "lfdr_tresh_viter"~  "S(X) == SC(FDR[v])",
                     Nom == "block" ~ "S(X) ==bgroup('{','1, ..., 200' ,'}')"), 
    Nom = factor(Nom, levels = c("S(X) == bgroup('{',p < 0.05,'}')","S(X) == SC(FDR[p])",
                                 "S(X) == SC(0.05)", "S(X) == bgroup('{',p < th,'}')",
                                 "S(X) == bgroup('{',viterbi == 1,'}')", 
                                 "S(X) == bgroup('{',viterbi[or] == 1,'}')",
                                 "S(X) == SC(FDR[v])",
                                 "S(X) == bgroup('{','viterbi = 1 & size > 4','}')", 
                                 "S(X) ==bgroup('{','1, ..., 200' ,'}')",
                                 "S(X) == H1",
                                 "S(X) ==H0~intersect()~ bgroup('{',p < 0.05 ,'}')", 
                                 "S(X) ==H0~intersect()~bgroup('{',viterbi == 1 ,'}')"))
    ) 
  return(Res_all) 
}



#' Title
#'
#' @param Result_hmm_boot 
#'
#' @return
#' @export
#'
#' @examples
res_FDR <- function(Result_hmm_boot){
  alpha <- Result_hmm_boot$al %>% unique()
  Res_all <-   Result_hmm_boot %>% 
    filter(Nom != "pval_tresh_min_size" & Size > 0) %>% 
    mutate(
    q_aldemi_high = map_dbl(Est_HMM_boot, function(x){
      mutate(x, eps =   (V_HMM_oracle_boot_aldemi - V_HMM_est_boot_aldemi)/Size_boot ) %>% 
        pull(eps) %>% 
        quantile(probs=(1 - alpha / 2), na.rm =TRUE, type = 3)
    }),
    q_aldemi_low = map_dbl(Est_HMM_boot, function(x){
      vec <- x %>% filter(Size_boot >0) %>%
        mutate( eps =   (V_HMM_small_oracle_boot_aldemi - V_HMM_small_est_boot_aldemi)/Size_boot) %>% 
        pull(eps) %>% sort()
      
      vec[max(floor(length(vec) * alpha/2),1)]
      # quantile(probs= (alpha / 2), na.rm =TRUE, type = 3)
    }),
    q_samesel = map2_dbl(Est_HMM_boot, V_HMM_est_aldemi, function(x,y){
      mutate(x, eps =  y -  V_HMM_est_boot_aldemi_samesel ) %>% 
        pull(eps) %>% 
        quantile(probs= (1- alpha / 2), na.rm =TRUE, type = 3)
    }),
    q_samesel_small = map2_dbl(Est_HMM_boot, V_HMM_small_est_aldemi, function(x,y){
      vec<-x %>% filter(Size_boot >0) %>%
        mutate(  eps =  y-  V_HMM_small_est_boot_aldemi_samesel ) %>% 
        pull(eps) %>% sort() 
      
      vec[max(floor(length(vec) * alpha/2),1)]
      # quantile(probs= (alpha / 2), na.rm =TRUE, type = 3)
    }),
    V_HMM_boot_samesel = V_HMM_est_aldemi + q_samesel,
    V_HMM_boot_samesel= V_HMM_boot_samesel /Size,
    V_HMM_small_boot_samesel = V_HMM_small_est_aldemi + q_samesel_small,
    V_HMM_small_boot_samesel= V_HMM_small_boot_samesel/Size,
    V_HMM_boot_compReal = FDR_est/Size + High_boot,
    V_HMM_small_boot_compReal = FDR_est/Size + Low_boot, 
    V_HMM_boot_qdemi = V_HMM_est_aldemi/Size + q_aldemi_high,
    V_HMM_small_boot_qdemi = V_HMM_small_est_aldemi/Size + q_aldemi_low, 
    V_HMM_boot_naif = map_dbl(Est_HMM_boot, ~pull(. , Real_boot) %>% 
                                quantile(probs = (1 - alpha), type = 3, na.rm =TRUE)),
    V_HMM_boot_naif= V_HMM_boot_naif/Size,
    V_HMM_small_boot_naif = map_dbl(Est_HMM_boot, function(x){
      vec <-  pull(x , Real_boot) %>% sort()
      vec[floor(length(vec) * alpha)]
    }),
    # quantile(probs = (alpha), type = 3, na.rm =TRUE)),
    V_HMM_small_boot_naif = V_HMM_small_boot_naif/Size,
    FDR_boot_naif = map_dbl(Est_HMM_boot, ~ pull(. , Real_boot) %>% 
                              mean()),
    V_simes = V_simes / Size, 
    FDR_boot_naif = FDR_boot_naif /Size, 
     FDR_est = FDR_est/Size,
    Nom = case_when( Nom == "pval_tresh" ~  "S(X) == bgroup('{',p < 0.05,'}')", 
                     Nom == "sel_viter_est" ~ "S(X) == bgroup('{',viterbi == 1,'}')",
                     Nom == "sel_viter_or" ~ "S(X) == bgroup('{',viterbi[or] == 1,'}')",
                     Nom == "sel_viter_min_size" ~ "S(X) == bgroup('{','viterbi = 1 & size > 4','}')",
                     Nom == "H0_and_pval_tresh_small" ~  "S(X) ==H0~intersect()~ bgroup('{',p < 0.05 ,'}')", 
                     Nom == "H0_and_sel_viter_est" ~ "S(X) ==H0~intersect()~bgroup('{',viterbi == 1 ,'}')",
                     Nom == "H1" ~  "S(X) == H1", 
                     Nom == "pval_ord" ~  "S(X) == bgroup('{',p < th,'}')",
                     Nom == "lfdr_tresh" ~  "S(X) == SC(0.05)",
                     Nom == "lfdr_tresh_pval" ~  "S(X) == SC(FDR[p])",
                     Nom == "lfdr_tresh_viter"~  "S(X) == SC(FDR[v])",
                     Nom == "block" ~ "S(X) ==bgroup('{','1, ..., 200' ,'}')"), 
    Nom = factor(Nom, levels = c("S(X) == bgroup('{',p < 0.05,'}')","S(X) == SC(FDR[p])",
                                 "S(X) == SC(0.05)", "S(X) == bgroup('{',p < th,'}')",
                                 "S(X) == bgroup('{',viterbi == 1,'}')", 
                                 "S(X) == bgroup('{',viterbi[or] == 1,'}')",
                                 "S(X) == SC(FDR[v])",
                                 "S(X) == bgroup('{','viterbi = 1 & size > 4','}')", 
                                 "S(X) ==bgroup('{','1, ..., 200' ,'}')",
                                 "S(X) == H1",
                                 "S(X) ==H0~intersect()~ bgroup('{',p < 0.05 ,'}')", 
                                 "S(X) ==H0~intersect()~bgroup('{',viterbi == 1 ,'}')"))
    ) 
  return(Res_all) 
}
