#' Title
#'
#' @param m 
#' @param A 
#' @param Pi 
#' @param f0x 
#' @param f1x 
#' @param x 
#' @param eps 
#' @param maxit 
#' @param h 
#'
#' @return
#' @export
#'
#' @examples
#' m <- 200
#' A <- matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T)
#' Pi <- c(0.8,0.2)
#' theta <- sim_markov(m, Pi , A)
#' x <- rep(0, m)
#' x[theta == 0] <- rnorm(sum(theta ==0))
#' x[theta == 1] <- rnorm(sum(theta ==1), 2, 1)
#' f0x <- dnorm(x)
#' f1x <- dnorm(x, 2,1)
#' Em_tot_01_approx(m, A, Pi, f0x, f1x,
#'  x, eps= 0.0001, maxit=1000, h= 0.3)
Em_tot_01_approx<- function(m, A, Pi, f0x, f1x,
   x, eps, maxit, h) {
 diff <- list(eps +1, eps +1,eps +1,eps +1)
   i <- 0
  y <- list()
  while(((max(unlist(diff)) > eps) & (i < maxit)))
  {
    # A_old = A
    # Pi_old = Pi
    # f0x_old = f0x
    # f1x_old = f1x
    # b_f = for_back(m, A, f0x, f1x, Pi)
    # sum_gamma = colSums(b_f$gamma[-m,])
    # sum_ksi = colSums(b_f$ksi[-m,])
    # Pi = b_f$gamma[1,]
    # A <- matrix(sum_ksi / rep(sum_gamma,2), byrow = FALSE, ncol = 2)
    # d0 <- density(x, weights = b_f$gamma[,1]/sum(b_f$gamma[,1]))
    # d1 <- density(x, weights = b_f$gamma[,2]/sum(b_f$gamma[,2]))
    # f1x <- approx(d1$x,d1$y,x)$y
    # f0x <- approx(d0$x,d0$y,x)$y
    # 
    # diff[1] = max(abs(A - A_old))
    # diff[2] = max(abs(Pi - Pi_old))
    # diff[3] = max(abs(f1x - f1x_old))
    # diff[4] = max(abs(f0x - f0x_old))
    y[[1]] <- y[[6]]
    y[[2]] <- y[[5]]$gamma[1,]
    y[[3]] <- y[[10]]
    y[[4]] <- y[[9]]
    y[[5]] <- for_back(m, y[[6]], f0x=y[[3]], f1x = y[[4]], y[[5]]$gamma[1,])
    y[[6]] <- matrix(colSums(y[[5]]$ksi[-m,]) / rep(colSums(y[[5]]$gamma[-m,]),2), byrow = FALSE, ncol = 2)
    y[[7]] <- density(x, weights = y[[5]]$gamma[,1]/sum(y[[5]]$gamma[,1]))
    y[[8]] <- density(x, weights = y[[5]]$gamma[,2]/sum(y[[5]]$gamma[,2]))
    y[[9]] <- approx(y[[8]]$x,y[[8]]$y,x)$y
    y[[10]] <- approx(y[[7]]$x,y[[7]]$y,x)$y
    
    diff[[1]] = max(abs( y[[6]] - y[[1]]))
    diff[[2]] = max(abs(y[[5]]$gamma[1,] - y[[2]]))
    diff[[3]] = max(abs(y[[9]] - y[[4]]))
    diff[[4]] = max(abs(y[[10]] - y[[3]]))
   i <- i+1
   gc()
  }
  
  b_f = for_back(m, A, f0x, f1x, Pi)
  return (list(A= A,Pi = Pi,
                            fw_bc_EM= b_f,
                            f1x = f1x,
                           f0x = f0x,
                            i = i))
}


#' Title
#'
#' @param m 
#' @param A 
#' @param Pi 
#' @param f0x 
#' @param f1x 
#' @param x 
#' @param eps 
#' @param maxit 
#' @param h 
#'
#' @return
#' @export
#'
#' @examples
#' m <- 200
#' A <- matrix(c(0.95, 0.05, 0.2, 0.80), 2, 2, byrow = T)
#' Pi <- c(0.8,0.2)
#' theta <- sim_markov(m, Pi , A)
#' x <- rep(0, m)
#' x[theta == 0] <- rnorm(sum(theta ==0))
#' x[theta == 1] <- rnorm(sum(theta ==1), 2, 1)
#' f0x <- dnorm(x)
#' f1x <- dnorm(x, 2,1)
#' Em_tot_approx(m, A, Pi, f0x, f1x,
#'  x, eps= 0.0001, maxit=1000, h= 0.3)
Em_tot_approx<- function(m, A, Pi, f0x, f1x,
                            x, eps, maxit, h) {
  diff <- rep(eps +1, 3)
  i <- 0
  y <- as.list(rep(0,8))
  y[[6]] <- A
  y[[3]] <- Pi
  y[[8]] <- f1x
  while(((max(unlist(diff)) > eps) & (i < maxit)))
  {
    # A_old = A
    # Pi_old = Pi
    # f1x_old = f1x
    # b_f = for_back(m, A, f0x, f1x, Pi)
    # sum_gamma = colSums(b_f$gamma[-m,])
    # sum_ksi = colSums(b_f$ksi[-m,])
    # Pi = b_f$gamma[1,]
    # A <- matrix(sum_ksi / rep(sum_gamma,2), byrow = FALSE, ncol = 2)
    # d0 <- density(x, weights = b_f$gamma[,1]/sum(b_f$gamma[,1]))
    # d1 <- density(x, weights = b_f$gamma[,2]/sum(b_f$gamma[,2]))
    # f1x <- approx(d1$x,d1$y,x)$y
    # 
    # diff[1] = max(abs(A - A_old))
    # diff[2] = max(abs(Pi - Pi_old))
    # diff[3] = max(abs(f1x - f1x_old))
    
    y[[1]] <- y[[6]]
    y[[2]] <-  y[[3]]
    y[[4]] <- y[[8]] 
    y[[5]] <- for_back(m,y[[1]], f0x, y[[4]], y[[3]])
    y[[3]] <- y[[5]]$gamma[1,]
    y[[6]] <- matrix(colSums(y[[5]]$ksi[-m,]) / rep(colSums(y[[5]]$gamma[-m,]),2), byrow = FALSE, ncol = 2)
    y[[7]] <- density(x, weights = y[[5]]$gamma[,2]/sum(y[[5]]$gamma[,2]))
    y[[8]] <- approx(y[[7]]$x,y[[7]]$y,x)$y
    
    diff[[1]] = max(abs(y[[6]] - y[[1]]))
    diff[[2]] = max(abs(y[[3]] - y[[2]]))
    diff[[3]] = max(abs(y[[8]] - y[[4]]))
    i <- i+1
    
    gc()
    
  }
  
  b_f = for_back(m, A, f0x, f1x, Pi)
  return (list(A= A,Pi = Pi,
               fw_bc_EM= b_f,
               f1x = f1x,
               f0x = f0x,
               i = i))
}