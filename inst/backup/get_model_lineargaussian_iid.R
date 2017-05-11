#'@rdname get_model_lineargaussian_iid
#'@title get_model_lineargaussian_iid
#'@description This implements a degenerate univariate linear Gaussian model with iid states
#'@export
get_model_lineargaussian_iid <- function(){
  model.lineargaussian = list()
  model.lineargaussian$df = 20
  model.lineargaussian$sigmaW2 = 0

  # dimension of parameter
  model.lineargaussian$dimtheta = 1
  model.lineargaussian$dimY = 1
  model.lineargaussian$dimX = 1

  # sampler from the prior distribution on parameters
  model.lineargaussian$rprior = function(Ntheta){
    nu = model.lineargaussian$df
    return (cbind(rinvchisq(Ntheta,nu,1)))
  }

  # prior distribution density on parameters
  model.lineargaussian$dprior = function(theta, log = TRUE){
    nu = model.lineargaussian$df
    sigmaV2 = theta[1]
    return (dinvchisq(sigmaV2,nu,1,log))
  }

  # sampler from the initial distribution of the states
  model.lineargaussian$rinitial = function(theta,N){
    sigmaW2 = model.lineargaussian$sigmaW2
    return (matrix(rnorm(N, mean = 0, sd = sqrt(sigmaW2)), nrow = N))
  }

  # sampler from the transition density of the states
  model.lineargaussian$rtransition = function(Xt,t,theta){
    sigmaW2 = model.lineargaussian$sigmaW2
    N = nrow(Xt)
    return (matrix(rnorm(N, mean = 0, sd = sqrt(sigmaW2)), nrow = N))
  }

  # density of the observations
  model.lineargaussian$dobs = function(Yt,Xt,t,theta,log = TRUE){
    sigmaV2 = theta[1]
    return (dnorm(Yt,mean = Xt,sd = sqrt(sigmaV2), log))
  }

  # first and second partial derivatives (k-th coordinate) of the observation log-density
  model.lineargaussian$derivativelogdobs = function(Yt,Xt,t,theta,k){
    sigmaV2 = theta[1]
    N = nrow(Xt)
    d1 = (Xt-matrix(Yt,nrow = N))/sigmaV2
    d2 = matrix(-1/sigmaV2,nrow = N)
    return (list(d1log = d1, d2log = d2))
  }

  # OPTIONAL: simulate observations
  model.lineargaussian$robs = function(Xt,t,theta){
    sigmaV2 = theta[1]
    N = nrow(Xt)
    return (Xt + rnorm(N, mean = 0, sd = sqrt(sigmaV2)))
  }

  return(model.lineargaussian)
}
