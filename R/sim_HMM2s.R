#' Title
#'
#' @param m
#' @param Pi
#' @param A
#'
#' @return
#' @export
#'
#' @examples
sim_markov <-function(m, Pi, A){

  theta <- rep(0, m)
  x <- rep(0, m)

  ## generating theta :
  theta[1] <- rbinom(1, 1, Pi[2])
  for (i in 2:m)
  { theta[i] <- (1-theta[i-1])*rbinom(1, 1, A[1, 2]) + theta[i-1]*rbinom(1, 1, A[2, 2])
  }

  return(theta)
}

#' Title
#'
#' @param m
#' @param Pi
#' @param A
#' @param f0
#' @param f1
#'
#' @return
#' @export
#'
#' @examples
sim_hmm_2states <-function(m, Pi, A, f0=c(0,1), f1)
{
  ## ARGUMENTS :
  # m: sample size
  # pi=(pi[1], pi[2]): initial state distribution
  # A=(a00, a01 \\ a10 a11): transition matrix
  # f0: parameter set for the null distribution.
  # f1: parameter set for the non-null distribution

  ## VALUES : list with
  # x: continuous observed data
  # theta: binary unobserved states

  ## DETAILS :
  #  theta is a markov chain with two states (0,1) characterized by initial distribution Pi and transition matrix A
  # xi are independent given theta with distribution                   p(xi|theta_i=0)=dnorm_f0(xi) and  p(xi|theta_i=1)=dnorm_f1(xi)

  theta <- rep(0, m)
  x <- rep(0, m)

  ## generating theta :
  theta[1] <- rbinom(1, 1, Pi[2])
  for (i in 2:m)
  { theta[i] <- (1-theta[i-1])*rbinom(1, 1, A[1, 2]) + theta[i-1]*rbinom(1, 1, A[2, 2])
  }

  ## generating the observations  #eviter boucle sur m :
  #for (i in 1:m)
  # { x[i]<- (1-theta[i])*rnorm(1, mean=f0[1], sd=f0[2]) + theta[i]*
  #rnorm(1, mean=f1[1], sd=f1[2])
  # }
  for (ind in 0:1)
  { nb.ind <- sum(theta==ind)
  f=(1-ind)*f0 + ind *f1  #ou f= list(f0,f1) puis rnorm(nb.ind, mean=f[ind][1], sd=f[ind][2])
  if (nb.ind>0) x[theta==ind] <- rnorm(nb.ind, mean=f[1], sd=f[2])
  }

  data<-list(theta=theta, x=x)
  return (data)
}


################################################################
#  forward-backward procedure : inspired by SUN and CAI
################################################################

#' Title
#'
#' @param x
#' @param Pi
#' @param A
#' @param f0
#' @param f1
#'
#' @return
#' @import Matrix
#' @export
#'
#' @examples
bwfw.hmm.2states <-function(x, Pi, A, f0=c(0,1), f1)
{
  ## ARGUMENTS
  # x=(x_1, ..., x_m): the observed data
  # Pi=(Pi[0], Pi[1]): the initial distribution
  # A=(A[0,0], A[0,1]\\ A[1,0], A[1,1]): the transition matrix
  # f0: the parameters for null distribution (mean, sd)
  # f1: the parameters for the non-null distribution

  ## VALUES
  # alpha: rescaled forward variables
  #alpha_i(j)=P(x_1, ..., x_i, theta_i =j)
  # beta: rescaled backward variables
  #beta_i(j)=P(x_i+1, ..., x_m | theta_i =j)
  # gamma: probabilities of hidden states :
  #gamma_i(j)=P(theta_i=j |x)
  # transition: rescaled transition variables :
  #transition_i(j,k)=P(theta_i =j,theta_i+1 =k|x)
  ## DETAILS
  # using the forward-backward procedure (Baum et al.)
  # the underflow problem was fixed by using the rescaled forward and backward variables : see Stamp (2004) and Rabiner(1989)


  m<-length(x)
  f0x<-dnorm(x, f0[1], f0[2])
  f1x<-dnorm(x, f1[1], f1[2])

  # ---- alpha: rescaled forward variables -----
  alpha<-matrix(NA, nrow=m, ncol=2)
  # scaling variable c_0
  c0<-rep(0, m)

  # initialize for i=1
  alpha[1, 1]<-Pi[1]*f0x[1]
  alpha[1, 2]<-Pi[2]*f1x[1]
  # rescaling alpha
  c0[1]<-1/sum(alpha[1, ])
  alpha[1, ]<-c0[1]*alpha[1, ]

  for (i in 1:(m-1))
  {
    alpha[i+1, 1] <-f0x[i+1] * (alpha[i, 1]*A[1, 1]+alpha[i, 2]*A[2, 1])
    alpha[i+1, 2] <-f1x[i+1] * (alpha[i, 1]*A[1, 2]+alpha[i, 2]*A[2, 2])
    # rescaling alpha
    c0[i+1] <- 1/sum(alpha[i+1, ])
    alpha[i+1, ]<-c0[i+1]*alpha[i+1, ]
  }

  #for (j in 0:1)
  #{ fx_j= (1-j)*f0x[1]+j*f1x[1]
  #  Pi_j=(1-j)*Pi[1]+j*Pi[2]
  #  alpha[1, j+1]<-Pi_j*fx_j


  # ---- beta: rescaled backward variables -----
  beta<-matrix(NA, nrow=m, ncol=2)

  # initialize for i=m.
  beta[m, 1]<-1
  beta[m, 2]<-1
  # rescaling beta using the same scaling factors as alpha
  beta[m, ]<-c0[m]*beta[m,]

  for (i in (m-1):1)
  {
    beta[i,1]<-A[1,1]*f0x[i+1]*beta[i+1,1]+A[1,2]*f1x[i+1]*beta[i+1,2]
    beta[i,2]<-A[2,1]*f0x[i+1]*beta[i+1,1]+A[2,2]*f1x[i+1]*beta[i+1,2]
    # rescaling beta
    beta[i, ]<-c0[i]*beta[i, ]
  }


  # -----gamma: probabilities of hidden states----------------

  gamma <- matrix(NA, nrow=m, ncol=2)
  for (i in 1:m)
  {  q1<-alpha[i, 1]*beta[i, 1]
  q2<-alpha[i, 2]*beta[i, 2]
  gamma[i,] <- c(q1,q2)/(q1+q2)
  }

  # LSI for state=0  : P(theta_i=0 |x)
  lsi=gamma[,1]


  # ---------- rescaled transition variables-----------------
  transition <-array(NA, c(2, 2, (m-1)))

  for (i in 1:(m-1))
  {
    denom<-0
    for (j in 0:1)
    {  for (k in 0:1)
    { fx<-(1-k)*f0x[i+1]+k*f1x[i+1]  #utiliser aussi pour alpha, beta?
    denom<-denom+alpha[i, j+1]*A[j+1, k+1]*fx*beta[i+1, k+1]
    }}

    for (j in 0:1)
    {  for (k in 0:1)
    { fx<-(1-k)*f0x[i+1]+k*f1x[i+1]
    transition[j+1, k+1, i] <- alpha[i, j+1]*A[j+1, k+1]*fx*beta[i+1, k+1]/denom
    }}
  }

  #--------- return the results of the bwfw proc-----------------
  bwfw.var<-list(alpha=alpha, beta=beta, gamma=gamma, transition=transition, lsi=lsi)

  return(bwfw.var)
}
