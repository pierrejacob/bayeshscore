#'@rdname get_model_iid_gaussian_unknown_variance
#'@title get_model_iid_gaussian_unknown_variance
#'@description Univariate iid Gaussian observations with unknown variance (degenerate linear gaussian SSM)
#'@export
get_model_iid_gaussian_unknown_variance <- function(){
  model = list()
  model$observation_type = 'continuous'
  model$df = 20
  model$sigmaW2 = 0

  # dimension of parameter
  model$dimtheta = 1
  model$dimY = 1
  model$dimX = 1

  # sampler from the prior distribution on parameters
  model$rprior = function(Ntheta){
    return (rbind(rinvchisq(Ntheta,model$df,1)))
  }

  # prior distribution density on parameters
  model$dprior = function(theta, log = TRUE){
    sigmaV2 = theta[1]
    return (dinvchisq(sigmaV2,model$df,1,log))
  }

  # sampler from the initial distribution of the states
  model$rinitial = function(theta,N){
    return (matrix(rnorm(N, mean = 0, sd = sqrt(model$sigmaW2)), ncol = N))
  }

  # sampler from the transition density of the states
  model$rtransition = function(Xt,t,theta){
    N = ncol(Xt)
    return (matrix(rnorm(N, mean = 0, sd = sqrt(model$sigmaW2)), ncol = N))
  }

  # density of the observations
  model$dobs = function(Yt,Xt,t,theta,log = TRUE){
    sigmaV2 = theta[1]
    return (dnorm(Yt,mean = Xt,sd = sqrt(sigmaV2), log))
  }

  # first and second partial derivatives of the observation log-density
  # The function is vectorized with respect to the states Xt (dimX by Nx), so that it outputs:
  # >> the jacobian (Nx by dimY matrix: each row is the transpose of the corresponding gradients row-wise)
  # >> the Hessian diagonals (Nx by dimY matrix: each row is the diagonal coeffs of the corresponding Hessian)
  model$derivativelogdobs = function(Yt,Xt,t,theta,dimY){
    sigmaV2 = theta[1]
    N = ncol(Xt)
    d1 = matrix((Xt-repeat_column(N,Yt))/sigmaV2,nrow = N, ncol = dimY)
    d2 = matrix(-1/sigmaV2,nrow = N, ncol = dimY)
    return (list(jacobian = d1, hessiandiag = d2))
  }
  # OPTIONAL: one-step predicitve density of the observation at time t given all the past from 1 to (t-1)
  # This relies on some Kalman filter (passed as a byproduct)
  model$dpredictive = function(observations,t,theta,log = TRUE){
    sigmaV2 = theta[1]
    return (dnorm(observations[,t],mean = 0,sd = sqrt(sigmaV2), log))
  }

  # OPTIONAL: simulate observations
  model$robs = function(Xt,t,theta){
    sigmaV2 = theta[1]
    N = ncol(Xt)
    return (matrix(Xt + rnorm(N, mean = 0, sd = sqrt(sigmaV2)),ncol = N))
  }

  return(model)
}
