#'@rdname get_model_lineargaussian
#'@title get_model_lineargaussian
#'@description This implements the univariate linear Gaussian model
#'@export
get_model_lineargaussian <- function(){
  model.lineargaussian = list()

  # dimension of parameter
  model.lineargaussian$dimtheta = 4
  model.lineargaussian$dimY = 1
  model.lineargaussian$dimX = 1

  # sampler from the prior distribution on parameters
  model.lineargaussian$rprior = function(Ntheta){
    phi = runif(Ntheta,0.1,0.9)
    psi = runif(Ntheta,0.5,1.5)
    sigmaV2 = runif(Ntheta,0.1,2)
    sigmaW2 = runif(Ntheta,1,10)
    return (cbind(phi,psi,sigmaV2,sigmaW2))
  }

  # prior distribution density on parameters
  model.lineargaussian$dprior = function(theta, log = TRUE){
    phi = theta[1]
    psi = theta[2]
    sigmaV2 = theta[3]
    sigmaW2 = theta[4]
    if (log==TRUE){
      return (dunif(phi,0.1,0.9,log)+dunif(psi,0.5,1.5,log)+dunif(sigmaV2,0.1,2,log)+dunif(sigmaW2,1,10,log))
    }
    else{
      return (dunif(phi,0.1,0.9,log)*dunif(psi,0.5,1.5,log)*dunif(sigmaV2,0.1,2,log)*dunif(sigmaW2,1,10,log))
    }
  }

  # sampler from the initial distribution of the states
  model.lineargaussian$rinitial = function(theta,N,additional_parameters){
    phi = theta[1]
    sigmaW2 = theta[4]
    return (matrix(rnorm(N, mean = 0, sd = sqrt((sigmaW2)/(1-phi^2))), nrow = N))
  }

  # sampler from the transition density of the states
  model.lineargaussian$rtransition = function(Xt,t,theta,additional_parameters){
    phi = theta[1]
    sigmaW2 = theta[4]
    N = nrow(Xt)
    return (matrix(phi*Xt + rnorm(N, mean = 0, sd = sqrt(sigmaW2)), nrow = N))
  }

  # density of the observations
  model.lineargaussian$dobs = function(Yt,Xt,t,theta,log = TRUE){
    psi = theta[2]
    sigmaV2 = theta[3]
    return (dnorm(Yt,mean = psi*Xt,sd = sqrt(sigmaV2), log))
  }

  # first and second partial derivatives (k-th coordinate) of the observation log-density
  model.lineargaussian$derivativelogdobs = function(Yt,Xt,t,theta,k){
    phi = theta[1]
    psi = theta[2]
    sigmaV2 = theta[3]
    sigmaW2 = theta[4]
    N = nrow(Xt)
    d1 = (phi*Xt-matrix(Yt,nrow = N))/sigmaV2
    d2 = matrix(-1/sigmaV2,nrow = N)
    return (list(d1log = d1, d2log = d2))
  }

  # OPTIONAL: simulate observations
  model.lineargaussian$robs = function(Xt,t,theta){
    psi = theta[2]
    sigmaV2 = theta[3]
    N = nrow(Xt)
    return (psi*Xt + rnorm(N, mean = 0, sd = sqrt(sigmaV2)))
  }

  # OPTIONAL: fixed parameter for simulation
  model.lineargaussian$theta = c(0.8,1,4,1)

  # OPTIONAL: initial mean and variance for state (for linear gaussian model only)
  model.lineargaussian$initialmean = 0
  model.lineargaussian$initialvar = (model.lineargaussian$theta[4])/(1-model.lineargaussian$theta[1]^2)
  return(model.lineargaussian)
}
