#'@rdname get_model_iid_gaussian_unknown_mean
#'@title get_model_iid_gaussian_unknown_mean
#'@description Univariate iid Gaussian observations with unknown mean (degenerate linear gaussian SSM)
#'@export
get_model_iid_gaussian_unknown_mean <- function(muprior,sigma2prior){
  model = list()
  model$observation_type = 'continuous'

  # dimension of parameter
  model$dimtheta = 1
  model$dimY = 1
  model$dimX = 1

  # sampler from the prior distribution on parameters
  model$rprior = function(Ntheta){
    return (rbind(rnorm(Ntheta,muprior,sqrt(sigma2prior))))
  }

  # prior distribution density on parameters
  model$dprior = function(theta, log = TRUE){
    return (dnorm(theta,muprior,sqrt(sigma2prior),log))
  }

  # sampler from the initial distribution of the states
  model$rinitial = function(theta,N){
    return (matrix(theta, ncol = N))
  }

  # sampler from the transition density of the states
  model$rtransition = function(Xt,t,theta){
    N = ncol(Xt)
    return (matrix(theta, ncol = N))
  }

  # density of the observations
  model$dobs = function(Yt,Xt,t,theta,log = TRUE){
    return (dnorm(Yt, theta, 1, log))
  }

  # first and second partial derivatives of the observation log-density
  # The function is vectorized with respect to the states Xt (dimX by Nx), so that it outputs:
  # >> the jacobian (Nx by dimY matrix: each row is the transpose of the corresponding gradients row-wise)
  # >> the Hessian diagonals (Nx by dimY matrix: each row is the diagonal coeffs of the corresponding Hessian)
  model$derivativelogdobs = function(Yt,Xt,t,theta,dimY){
    N = ncol(Xt)
    d1 = t((matrix(theta,ncol=N)-repeat_column(N,Yt)))
    d2 = matrix(-1,nrow = N, ncol = dimY)
    return (list(jacobian = d1, hessiandiag = d2))
  }
  # OPTIONAL: one-step predicitve density of the observation at time t given all the past from 1 to (t-1)
  # This relies on some Kalman filter (passed as a byproduct)
  model$dpredictive = function(observations,t,theta,log = TRUE){
    return (dnorm(observations[,t], theta, 1, log))
  }

  # OPTIONAL: simulate observations
  model$robs = function(Xt,t,theta){
    N = ncol(Xt)
    return (matrix(rnorm(N, theta, 1),ncol = N))
  }

  return(model)
}
