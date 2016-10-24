#'@rdname get_model_pz4
#'@title get_model_pz4
#'@description This implements the PZ model
#'@export
get_model_pz4 <- function(){

  model.pz = list()

  # dimension of parameter, observations, and states
  model.pz$dimtheta = 4
  model.pz$dimY = 1
  model.pz$dimX = 2

  # sampler from the prior distribution on parameters
  model.pz$rprior = function(Ntheta){
    theta1 <- runif(Ntheta)
    theta2 <- runif(Ntheta)
    theta3 <- runif(Ntheta)
    theta4 <- runif(Ntheta)
    return (cbind(theta1,theta2,theta3,theta4))
  }

  # density the prior distribution on parameters
  model.pz$dprior = function(theta, log = TRUE){
    density <- dunif(theta[1],0,1, log = TRUE) + dunif(theta[2],0,1, log = TRUE) +
      dunif(theta[3],0,1, log = TRUE) + dunif(theta[4],0,1, log = TRUE)
    if (log){
      return(density)
    } else {
      return(exp(density))
    }
  }

  # sampler from the initial distribution of the states
  model.pz$rinitial = function(theta, N){
    return (matrix(rlnorm(N*2,meanlog = log(2), sdlog = 1), ncol = 2))
  }

  model.pz$rtransition = function(Xt, t, theta){
    N <- nrow(Xt)
    alphas <- theta[1] + theta[2] * rnorm(N, 0, 1)
    xparticles <- pz_transition(Xt, alphas, t-1, c(theta[3:4], 0.1, 0.1))
    return(xparticles)
  }

  # density of the observations
  model.pz$dobs = function(Yt,Xt,t,theta, log = TRUE){
    # tau = theta[3]
    # n = 1/tau
    # p = 1/(1+tau*Xt)
    # return (dnbinom(Yt[1],size = n,prob = p)*dnbinom(Yt[2],size = n,prob = p))
    return(dnorm(Yt[1], mean = log(Xt[,1]), sd = 0.2, log = log))
  }

  model.pz$derivativelogdobs = function(Yt,Xt,t,theta,k){
    N = nrow(Xt)
    d1 = (log(Xt[,1])-matrix(Yt, N))/(0.2^2)
    d2 = matrix(-1/(0.2^2),nrow = N)
    return (list(d1log = d1, d2log = d2))
  }
  model.pz$robs = function(Xt,t,theta){
    N = nrow(Xt)
    return (log(Xt[,1]) + rnorm(N, mean = 0, sd = 0.2))
  }
  return(model.pz)
}
