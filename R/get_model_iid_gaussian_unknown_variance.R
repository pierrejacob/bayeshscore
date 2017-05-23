#'@rdname get_model_iid_gaussian_unknown_variance
#'@title get_model_iid_gaussian_unknown_variance
#'@description Univariate iid Gaussian observations with unknown variance (degenerate linear gaussian SSM)
#'@export
get_model_iid_gaussian_unknown_variance <- function(priordf){
  model = list()
  model$observation_type = 'continuous'

  # dimension of parameter
  model$dimtheta = 1
  model$dimY = 1

  # sampler from the prior distribution on parameters
  model$rprior = function(Ntheta){
    return (rbind(rinvchisq(Ntheta,priordf,1)))
  }

  # prior distribution density on parameters
  model$dprior = function(theta, log = TRUE){
    return (dinvchisq(theta,priordf,1,log))
  }

  # one-step predicitve density of the observation at time t given all the past from 1 to (t-1)
  model$dpredictive = function(observations,t,theta,log = TRUE){
    return (dnorm(observations[,t],mean = 0,sd = sqrt(theta), log))
  }

  # OPTIONAL: simulate observations
  model$robs = function(Xt,t,theta){
    sigmaV2 = theta[1]
    N = ncol(Xt)
    return (matrix(Xt + rnorm(N, mean = 0, sd = sqrt(sigmaV2)),ncol = N))
  }

  return(model)
}
