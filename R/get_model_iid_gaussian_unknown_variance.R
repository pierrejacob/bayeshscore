#'@rdname get_model_iid_gaussian_unknown_variance
#'@title get_model_iid_gaussian_unknown_variance
#'@description Univariate iid Gaussian observations with unknown variance (degenerate linear gaussian SSM)
#'@export
get_model_iid_gaussian_unknown_variance <- function(nu0, sigma02){
  model = list()
  model$observation_type = 'continuous'

  # dimension of parameter
  model$dimtheta = 1
  model$dimY = 1
  # fix some known parameters
  model$mu = 0

  # sampler from the prior distribution on parameters
  model$rprior = function(Ntheta){
    return (rbind(rinvchisq(Ntheta,nu0,sigma02)))
  }

  # prior distribution density on parameters
  model$dprior = function(theta, log = TRUE){
    return (dinvchisq(theta,nu0,sigma02,log))
  }

  # one-step predicitve density of the observation at time t given all the past from 1 to (t-1)
  model$dpredictive = function(observations,t,theta,log = TRUE){
    return (dnorm(observations[,t],mean = model$mu,sd = sqrt(theta), log))
  }

  # OPTIONAL: derivatives of the predicitve density with respect to the observation at time t
  # inputs: observations (dimY by T matrix, with T >= t), time index t (int), theta (single vector),
  #         byproduct (OPTIONAL: auxiliary object needed to compute likelihood, e.g. Kalman filter)
  # outputs: list with the following fields
  # jacobian >> the transpose of the gradient (1 by dimY)
  # hessiandiag >> the Hessian diagonal coefficients (1 by dimY)
  # NB: if missing, this field is automatically filled with numerical derivatives
  # via set_default_model in util_default.R)
  model$derivativelogdpredictive = function(observations,t,theta,byproduct) {
    deriv1 <- -(observations[,t]-model$mu)/theta
    deriv2 <- -1/theta
    return (list(jacobian = matrix(deriv1, 1, 1), hessiandiag = matrix(deriv2, 1, 1)))
  }

  # OPTIONAL: simulate observations
  model$robs = function(nobservations,theta){
    return (matrix(rnorm(nobservations, mean = model$mu, sd = sqrt(theta)),ncol = nobservations))
  }

  return(model)
}
