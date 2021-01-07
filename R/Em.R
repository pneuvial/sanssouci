Em <- function(m, A ,
               Pi,  f0x, f1x,
               x, eps = 0.0001,
               maxit =1000, h = 0.3, known_f0){
  
  if(known_f0){
    Em_tot(m, A ,
              Pi,  f0x, f1x,
              x, eps,
              maxit, h)
  }else{
    Em_tot_01(m, A ,
              Pi,  f0x, f1x,
              x, eps,
              maxit, h)   
  }

}