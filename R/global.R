global <- function(x, mu_0, sd_0, h, K, seuil_lis, seuil_p,type_sel, al ){
  m <- nrow(df())
  pval <- 2 * (1 - pnorm(abs(x)))
  f0x <- dnorm(x, mu_0, sd_0)
  pi0_hat <- min(sum(pval > 0.8) / (m * 0.2), 0.97)
  f_hatx <- x %>%
    map_dbl( ~f_hatK(x, ., h = h,K) )
  f1x_est <-  f1x_hat(f0x, f_hatx, pi0_hat)
  f1x_est[f1x_est <= 0] <- min(f0x)
  
  a <- runif(1, 0.7, 0.98)
  b <- 1 - a
  c <- runif(1, 0.1, 0.3)
  d <- 1 - c
  
  Em <- Em_tot(m, A = matrix(c(a, b, c, d), byrow = TRUE, ncol=2),
               Pi= c(pi0_hat, 1 -pi0_hat),  f0x = f0x, 
               f1x = f1x_est,
               x , eps = 0.0001,
               maxit =1000, h = h)
  Em$Pis <-  lapply(2:m, function(i){
    get_A( m,alpha = Em$fw_bc_EM$alpha, beta = Em$fw_bc_EM$beta,
           Em$A, f0x,
           Em$f1x, i = i)
    
    Selection <-   tibble(Type_sel = type_sel) %>% 
      mutate(Sel = map(Type_sel,~selection(., Em, 
                                           list(x = x, m =m , pval = pval, f0x = f0x),
                                           seuil_lis, seuil_p) ), 
             Size = map_dbl(Sel, ~length(.)))
    
    quantiles <- Selection %>% mutate(quantiles =  
              map(Sel, ~get_quantiles(sel = .,
                                      li0 = Em$fw_bc_EM$gamma[,1],
                   Pis = Em$Pis, f0x = f0x, f1x = Em$f1x)
                ))

    
    
    V <- 
      quantiles %>%  
        mutate(V_est = map2_dbl(Sel, quantiles,~borne(
          type_borne = "HMM", sel = .x, a = .y, 
          alpha = al)
        ))
    
    
  })
}