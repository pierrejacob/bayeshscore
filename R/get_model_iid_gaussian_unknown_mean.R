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

  # fix some known parameters
  model$sigma2 = 1

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
    return (dnorm(observations[,t], theta, sqrt(model$sigma2), log))
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
    deriv1 <- -(observations[,t]-theta)/model$sigma2
    deriv2 <- -1/model$sigma2
    return (list(jacobian = matrix(deriv1, 1, 1), hessiandiag = matrix(deriv2, 1, 1)))
  }

  # OPTIONAL: simulate observations
  model$robs = function(nobservations,theta){
    return (matrix(rnorm(nobservations, theta, sqrt(model$sigma2)),ncol = nobservations))
  }

  # OPTIONAL: simulate Ny draws of y_t given theta and the past y_1 to y_(t-1)
  # (with the convention y_0 = NULL)
  # outputs: matrix of Ny draws of Yt given theta and past (dimY by Ny matrix)
  model$rpredictive = function(Ny,t,theta,y_past){
    return (matrix(rnorm(Ny, theta, sqrt(model$sigma2)),ncol = Ny))
  }

  return(model)
}
