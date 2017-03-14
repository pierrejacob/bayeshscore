#'@rdname get_model_simplerlineargaussian
#'@title get_model_simplerlineargaussian
#'@description This implements the univariate linear Gaussian model
#'@export
get_model_simplerlineargaussian <- function(){
  model.lineargaussian = list()

  # dimension of parameter
  model.lineargaussian$dimtheta = 2
  model.lineargaussian$dimY = 1
  model.lineargaussian$dimX = 1

  # sampler from the prior distribution on parameters
  model.lineargaussian$rprior = function(Ntheta){
    phi = runif(Ntheta,0.1,0.9)
    # psi = runif(Ntheta,0.5,1.5)
    # sigmaV2 = runif(Ntheta,0.1,10)
    sigmaW2 = runif(Ntheta,0.1,10)
    return (cbind(phi,sigmaW2))
  }

  # prior distribution density on parameters
  model.lineargaussian$dprior = function(theta, log = TRUE){
    phi = theta[1]
    sigmaW2 = theta[2]
    if (log==TRUE){
      return (dunif(phi,0.1,0.9,log)+dunif(sigmaW2,0.1,10,log))
    }
    else{
      return (dunif(phi,0.1,0.9,log)*dunif(sigmaW2,0.1,10,log))
    }
  }

  # sampler from the initial distribution of the states
  model.lineargaussian$rinitial = function(theta,N){
    phi = theta[1]
    sigmaW2 = theta[2]
    return (matrix(rnorm(N, mean = 0, sd = sqrt((sigmaW2)/(1-phi^2))), nrow = N))
  }

  # sampler from the transition density of the states
  model.lineargaussian$rtransition = function(Xt,t,theta){
    phi = theta[1]
    sigmaW2 = theta[2]
    N = nrow(Xt)
    return (matrix(phi*Xt + rnorm(N, mean = 0, sd = sqrt(sigmaW2)), nrow = N))
  }

  # density of the observations
  model.lineargaussian$dobs = function(Yt,Xt,t,theta,log = TRUE){
    psi = 1
    sigmaV2 = 0.8
    return (dnorm(Yt,mean = psi*Xt,sd = sqrt(sigmaV2), log))
  }

  # first and second partial derivatives (k-th coordinate) of the observation log-density
  model.lineargaussian$derivativelogdobs = function(Yt,Xt,t,theta,k){
    phi = theta[1]
    sigmaW2 = theta[2]
    psi = 1
    sigmaV2 = 0.8
    N = nrow(Xt)
    d1 = (phi*Xt-matrix(Yt,nrow = N))/sigmaV2
    d2 = matrix(-1/sigmaV2,nrow = N)
    return (list(d1log = d1, d2log = d2))
  }

  # OPTIONAL: simulate observations
  model.lineargaussian$robs = function(Xt,t,theta){
    psi = 1
    sigmaV2 = 0.8
    N = nrow(Xt)
    return (psi*Xt + rnorm(N, mean = 0, sd = sqrt(sigmaV2)))
  }

  # OPTIONAL: fixed parameter for simulation
  model.lineargaussian$theta = c(0.8,0.9)

  # OPTIONAL: initial mean and variance for state (for linear gaussian model only)
  model.lineargaussian$initialmean = 0
  model.lineargaussian$initialvar = (model.lineargaussian$theta[2])/(1-model.lineargaussian$theta[1]^2)
  return(model.lineargaussian)
}
