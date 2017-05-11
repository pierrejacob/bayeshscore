#'@rdname get_model_lineargaussian_discreteprior
#'@title get_model_lineargaussian_discreteprior
#'@description This implements a univariate linear Gaussian model with discrete prior
#'@export
get_model_lineargaussian_discreteprior <- function(){
  model.lineargaussian = list()
  model.lineargaussian$supportprior = seq(0.1,3,0.01)
  model.lineargaussian$phi = 0.5
  model.lineargaussian$psi = 0.5
  model.lineargaussian$sigmaW2 = 1

  # dimension of parameter
  model.lineargaussian$dimtheta = 1
  model.lineargaussian$dimY = 1
  model.lineargaussian$dimX = 1

  # sampler from the prior distribution on parameters
  model.lineargaussian$rprior = function(Ntheta){
    support = model.lineargaussian$supportprior
    sigmaV2 = sample(support,Ntheta,replace = TRUE)
    return (rbind(sigmaV2))
  }

  # prior distribution density on parameters
  model.lineargaussian$dprior = function(theta, log = TRUE){
    support = model.lineargaussian$supportprior
    sigmaV2 = theta[1]
    if (sigmaV2 %in% support) {
      n = length(support)
      if (log==TRUE){
        return (-log(n))
      }
      else {
        return (1/n)
      }
    }
    else {
      if (log==TRUE){
        return (-Inf)
      }
      else {
        return (0)
      }
    }
  }

  # sampler from the initial distribution of the states
  model.lineargaussian$rinitial = function(theta,N){
    phi = model.lineargaussian$phi
    sigmaW2 = model.lineargaussian$sigmaW2
    return (matrix(rnorm(N, mean = 0, sd = sqrt((sigmaW2)/(1-phi^2))), ncol = N))
  }

  # sampler from the transition density of the states
  model.lineargaussian$rtransition = function(Xt,t,theta){
    phi = model.lineargaussian$phi
    sigmaW2 = model.lineargaussian$sigmaW2
    N = ncol(Xt)
    return (matrix(phi*Xt + rnorm(N, mean = 0, sd = sqrt(sigmaW2)), ncol = N))
  }

  # density of the observations
  model.lineargaussian$dobs = function(Yt,Xt,t,theta,log = TRUE){
    psi = model.lineargaussian$psi
    sigmaV2 = theta[1]
    return (dnorm(Yt,mean = psi*Xt,sd = sqrt(sigmaV2), log))
  }

  # first and second partial derivatives (k-th coordinate) of the observation log-density
  model.lineargaussian$derivativelogdobs = function(Yt,Xt,t,theta,k){
    phi = model.lineargaussian$phi
    psi = model.lineargaussian$psi
    sigmaV2 = theta[1]
    sigmaW2 = model.lineargaussian$sigmaW2
    N = ncol(Xt)
    d1 = (phi*Xt-matrix(Yt,ncol = N))/sigmaV2
    d2 = matrix(-1/sigmaV2,ncol = N)
    return (list(d1log = d1, d2log = d2))
  }

  # OPTIONAL: simulate observations
  model.lineargaussian$robs = function(Xt,t,theta){
    psi = model.lineargaussian$psi
    sigmaV2 = theta[1]
    N = ncol(Xt)
    return (psi*Xt + rnorm(N, mean = 0, sd = sqrt(sigmaV2)))
  }

  # OPTIONAL: initial mean and variance for state (for linear gaussian model only)
  model.lineargaussian$initialmean = 0
  model.lineargaussian$initialvar = (model.lineargaussian$sigmaW2)/(1-model.lineargaussian$phi^2)
  return(model.lineargaussian)
}
