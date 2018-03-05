##################################################################################################
# Generic template to define a model
##################################################################################################
model = list()
# Type of observations (string): "continuous" or "discrete"
model$observation_type = ...
# Dimension of parameter, observations, and possibly latent states (int)
model$dimtheta = ...
model$dimY = ...
model$dimX = ...
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# If the observations are discrete, lower and upper bounds for each component of the observations
# are required to compute the Hyvarinen score.
# model$lower and model$upper are vectors of length model$dimY, such that model$lower[j] and
# model$upper[j] are respectively the lower and upper bound defining the support of the j-th
# component of the observations
# Note: if no bounds are provided, the lower/upper bounds are set to -Inf/Inf by default
model$lower = rep(-Inf, model$dimY)
model$upper = rep(Inf, model$dimY)
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------

# Sampler from the prior distribution on parameters
# inputs: Ntheta (int)
# outputs: matrix (dimtheta by Ntheta) of prior draws
model$rprior = function(Ntheta){
  ...
  return (...)
}
# prior density on parameters
# inputs: theta (single vector of length dimtheta), log (TRUE by default)
# outputs: prior (log)-density theta (double)
model$dprior = function(theta, log = TRUE){
  lp = ...
  if (log==TRUE) {return (lp)}
  else {return (exp(lp))}
}
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Note: to use SMC, one may specify either the likelihood OR the one-step ahead predictive
# (one is automatically filled given the other, via set_default_model in util_default.R)
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# OPTIONAL: likelihood of the observations from 1 to t given theta
# inputs: observations (dimY by T matrix, with T >= t), time index t (int), theta (single vector of
#         length dimtheta), byproduct (OPTIONAL: auxiliary object needed to compute likelihood, e.g. Kalman filter),
#         log (TRUE by default)
# outputs: log-likelihood (or likelihood) of the observations from time 1 to t given theta (double)
# WARNING: must be an explicit function of the observation at time t to allow the
# computation of the derivative of the log-predictive density
model$likelihood = function(observations,t,theta,byproduct,log = TRUE){
  ll = ...
  if (log) {return(ll)}
  else {return(exp(ll))}
}
# OPTIONAL: one-step predicitve density of the observation at t given the past from 1 to (t-1) and theta
# inputs: observations (dimY by T matrix, with T >= t), time index t (int), theta (single vector of length dimtheta),
#         byproduct (OPTIONAL: auxiliary object needed to compute likelihood, e.g. Kalman filter),
#         log (TRUE by default)
# outputs: log-likelihood of the observations from time 1 to t given theta (double)
# WARNING: must be an explicit function of the observation at time t to allow the
# computation of the derivative of the log-predictive density
model$dpredictive = function(observations,t,theta,byproduct,log = TRUE){
  lp = ...
  if (log) {return(lp)}
  else {return(exp(lp))}
}


# OPTIONAL: Sampler from the one-step predictive given theta
# (i.e. y_t given theta and past y_1 to y_(t-1), with the convention y_0 = NULL)
# inputs: number of draws Ny, time t, parameter theta (single vector of length dimtheta),
# past Y (dimY by (t-1) matrix, or NULL if t=1)
# outputs: matrix of Ny draws of Yt given theta and past (dimY by Ny matrix)
# The following sampler is only relevant if one plans to simulate data from the model,
# or use kernel density estimator to compute the Hyvarinen score
model$rpredictive = function(Ny,t,theta,y_past){
  ...
  return (...)
}

# OPTIONAL: derivatives of the predicitve density with respect to the observation at time t
# inputs: observations (dimY by T matrix, with T >= t), time index t (int), theta (single vector of length dimtheta),
#         byproduct (OPTIONAL: auxiliary object needed to compute likelihood, e.g. Kalman filter)
# outputs: list with the following fields
# jacobian >> the transpose of the gradient (1 by dimY)
# hessiandiag >> the Hessian diagonal coefficients (1 by dimY)
# NB: if missing, this field is automatically filled with numerical derivatives
# via set_default_model in util_default.R)
model$derivativelogdpredictive = function(observations,t,theta,byproduct) {
  ...
  return (list(jacobian = ..., hessiandiag = ...))
}
# OPTIONAL: initialize byproducts if needed (e.g. Kalman filters, etc ...)
# inputs: theta (single vector of length dimtheta), observations (dimY by T matrix)
# outputs: byproduct associated with theta
model$initialize_byproducts = function(theta, observations){
  ...
  return(...)
}
# OPTIONAL: update byproducts if needed (e.g. Kalman filters, etc ...)
# inputs: current_byproduct at time t, time t, theta (single vector of length dimtheta),
# observations (dimY by T matrix)
# outputs: updated byproduct
model$update_byproduct = function(current_byproduct, t, theta, observations){
  ...
  return(...)
}
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Note: if no likelihood nor predictive is provided, the method will be SMC2, which requires
# specifying the transition kernel and the observation density
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Sampler from the initial distribution of the latent states
# inputs: theta (single vector of length dimtheta), Nx (int)
# outputs: matrix (dimX by Nx) of latent states
model$rinitial = function(theta,Nx){
  ...
  return (...)
}
# Sampler from the transition distribution of the latent states
# inputs: current states Xs at time (t-1) (dimX by Nx matrix), time t (int), theta (single vector of length dimtheta)
# outputs: updated states (dimX by Nx)
model$rtransition = function(Xs,t,theta){
  Nx = ncol(Xs)
  ...
  return (...)
}
# observation density
# inputs: single observation Yt (dimY by 1), states Xts (dimX by Nx), time t,
# theta (single vector of length dimtheta), log (TRUE by default)
# outputs: observation (log)-densities ("vectorized" with respect to the states Xt)
model$dobs = function(Yt,Xts,t,theta,log = TRUE){
  ...
  return (...)
}
# OPTIONAL: first and second partial derivatives of the observation log-density
# The function is vectorized with respect to the states Xts (dimX by Nx)
# inputs: single observation Yt (dimY by 1), states Xts (dimX by Nx), time t, theta (single vector of length dimtheta)
# outputs: list with the following fields
# >> the jacobian (Nx by dimY matrix: each row is the transpose of the corresponding gradients row-wise)
# >> the Hessian diagonals (Nx by dimY matrix: each row is the diagonal coeffs of the corresponding Hessian)
# NB: if missing, this field is automatically filled with numerical derivatives
# via set_default_model in util_default.R)
model$derivativelogdobs = function(Yt,Xts,t,theta){
  ...
  return (list(jacobian = ..., hessiandiag = ..))
}
# OPTIONAL: sampler from the observation distribution
# inputs: single state Xt (dimX by 1), time t, theta (single vector of length dimtheta)
# outputs: single observation (dimY by 1 matrix)
model$robs = function(Xt,t,theta){
  ...
  return (...)
}
