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

  # sampler from the prior distribution on parameters
  model$rprior = function(Ntheta){
    return (rbind(rnorm(Ntheta,muprior,sqrt(sigma2prior))))
  }

  # prior distribution density on parameters
  model$dprior = function(theta, log = TRUE){
    return (dnorm(theta,muprior,sqrt(sigma2prior),log))
  }

  # one-step predicitve density of the observation at time t given all the past from 1 to (t-1)
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
