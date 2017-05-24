#'@rdname get_model_iid_lognormal
#'@title get_model_iid_lognormal
#'@description Univariate iid lognormal observations
#'@export
get_model_iid_lognormal <- function(mu0,sigma02,a,b){
  model = list()
  # Type of observations (string): "continuous" or "discrete"
  model$observation_type = "continuous"
  # Dimension of parameter, observations, and possibly latent states (int)
  model$dimtheta = 2
  model$dimY = 1
  # Sampler from the prior distribution on parameters
  # inputs: Ntheta (int)
  # outputs: matrix (dimtheta by Ntheta) of prior draws
  model$rprior = function(Ntheta){
    mus <- rnorm(Ntheta, mu0, sqrt(sigma02))
    sigmas_square <- rinvgamma(Ntheta, a, b)
    return(rbind(mus, sigmas_square))
  }
  # prior density on parameters
  # inputs: theta (single vector), log (TRUE by default)
  # outputs: prior (log)-density theta (double)
  model$dprior = function(theta, log = TRUE){
    lp = dnorm(theta[1],mu0,sqrt(sigma02),TRUE) + dinvgamma(theta[2],a,b,TRUE)
    if (log==TRUE) {return (lp)}
    else {return (exp(lp))}
  }
  #----------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------
  # Note: to use SMC, one may specify either the likelihood OR the one-step ahead predictive
  # (one is automatically filled given the other, via set_default_model in util_default.R)
  #----------------------------------------------------------------------------------------------------
  #----------------------------------------------------------------------------------------------------
  # OPTIONAL: one-step predicitve density of the observation at t given the past from 1 to (t-1) and theta
  # inputs: observations (dimY by T matrix, with T >= t), time index t (int), theta (single vector),
  #         byproduct (OPTIONAL: auxiliary object needed to compute likelihood, e.g. Kalman filter),
  #         log (TRUE by default)
  # outputs: log-likelihood of the observations from time 1 to t given theta (double)
  # WARNING: must be an explicit function of the observation at time t to allow the
  # computation of the derivative of the log-predictive density
  model$dpredictive = function(observations,t,theta,byproduct,log = TRUE){
    y <- observations[,t]
    lp = -0.5*log(2*pi*theta[2]) - log(y) - 0.5*(log(y) - theta[1])^2 / theta[2]
    if (log) {return(lp)}
    else {return(exp(lp))}
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
    y <- observations[,t]
    deriv1 <- (-1/y) * (1 + (log(y) - theta[1]) / theta[2])
    deriv2 <- (1/(y^2)) * (1 + (log(y) - theta[1] - 1) / theta[2])
    return (list(jacobian = matrix(deriv1, 1, 1), hessiandiag = matrix(deriv2, 1, 1)))
  }
  return(model)
}
