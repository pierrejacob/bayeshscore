#------------------------------------------------------------------------------#
#-------------- Some useful functions to fill in missing arguments  -----------#
#------------------------------------------------------------------------------#
#'@rdname set_default_algorithmic_parameters
#'@title set_default_algorithmic_parameters
#'@description Set the missing algorithmic parameters to their default values
#'@export
set_default_algorithmic_parameters = function(algorithmic_parameters){
  if (is.null(algorithmic_parameters$adaptNx)) {algorithmic_parameters$adaptNx = TRUE}
  if (is.null(algorithmic_parameters$min_acceptance_rate)) {algorithmic_parameters$min_acceptance_rate = 0.20}
  if (is.null(algorithmic_parameters$ess_threshold)) {algorithmic_parameters$ess_threshold = 0.5}
  if (is.null(algorithmic_parameters$nmoves)) {algorithmic_parameters$nmoves = 1}
  if (is.null(algorithmic_parameters$store)) {algorithmic_parameters$store = FALSE}
  if (is.null(algorithmic_parameters$progress)) {algorithmic_parameters$progress = FALSE}
  if (is.null(algorithmic_parameters$verbose)) {algorithmic_parameters$verbose = FALSE}
  if (is.null(algorithmic_parameters$resampling)) {
    algorithmic_parameters$resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1))
  }
  return(algorithmic_parameters)
}
#------------------------------------------------------------------------------#
#'@rdname set_default_model
#'@title set_default_model
#'@description Set the missing field of a model to their default values
#'@export
set_default_model = function(model){
  # Define the derivatives of the observation log-density via numerical derivation (cf. numDeriv)
  # The function is vectorized with respect to the states Xt (dimX by Nx), so that it outputs:
  # >> the jacobian (Nx by dimY matrix: each row is the transpose of the corresponding gradients row-wise)
  # >> the Hessian diagonals (Nx by dimY matrix: each row is the diagonal coeffs of the corresponding Hessian)
  if (is.null(model$derivativelogdobs)){
    model$derivativelogdobs = function(Yt,Xt,t,theta,dimY){
      logdobs = function(y){model$dobs(y,Xt,t,theta,log = TRUE)}
      derivatives = get_derivatives_from_genD(genD(logdobs,Yt)$D,dimY)
      return (list(jacobian = derivatives$jacobian, hessiandiag = derivatives$hessiandiag))
    }
  }
  # Define the one-step predicitve density of the observation at time t given all the past from time 1 to (t-1)
  # This is only relevant when the likelihood is available
  if (!is.null(model$likelihood)){
    if (is.null(model$dpredictive)){
      # Some models require keeping track of some byproducts (e.g. Kalman filter recursion)
      if (is.null(model$initialize_byproducts)){
        model$dpredictive = function(observations,t,theta,log = TRUE){
          if (t==1){
            return(model$likelihood(observations,t,theta,log))
          } else {
            if (log){
              return(model$likelihood(observations,t,theta,log)-model$likelihood(observations,t-1,theta,log))
            } else {
              return(model$likelihood(observations,t,theta,log)/model$likelihood(observations,t-1,theta,log))
            }
          }
        }
      } else {
        model$dpredictive = function(observations,t,theta,byproduct,log = TRUE){
          if (t==1){
            return(model$likelihood(observations,t,theta,byproduct,log))
          } else {
            if (log){
              return(model$likelihood(observations,t,theta,byproduct,log)-model$likelihood(observations,t-1,theta,byproduct,log))
            } else {
              return(model$likelihood(observations,t,theta,byproduct,log)/model$likelihood(observations,t-1,theta,byproduct,log))
            }
          }
        }
      }
    }
  }
  # Define the derivatives of the predictive log-density via numerical derivation (cf. numDeriv)
  # The function outputs:
  # >> the transpose of the gradient (1 by dimY)
  # >> the Hessian diagonal coefficients (1 by dimY)
  if (is.null(model$derivativelogdpredictive)){
    if (is.null(model$initialize_byproducts)){
      model$derivativelogdpredictive = function(observations,t,theta,dimY){
        if (t==1){
          logpredictive = function(y) {model$dpredictive(y,t,theta,log = TRUE)}
        } else {
          logpredictive = function(y) {model$dpredictive(cbind(observations[1:(t-1)],y),t,theta,log = TRUE)}
        }
        derivatives = get_derivatives_from_genD(genD(logdobs,observations[,t])$D,dimY)
        return (list(jacobian = derivatives$jacobian, hessiandiag = derivatives$hessiandiag))
      }
    } else {
      model$derivativelogdpredictive = function(observations,t,theta,byproduct,dimY){
        if (t==1){
          logpredictive = function(y) {model$dpredictive(y,t,theta,byproduct,log = TRUE)}
        } else {
          logpredictive = function(y) {model$dpredictive(cbind(observations[1:(t-1)],y),t,theta,byproduct,log = TRUE)}
        }
        derivatives = get_derivatives_from_genD(genD(logdobs,observations[,t])$D,dimY)
        return (list(jacobian = derivatives$jacobian, hessiandiag = derivatives$hessiandiag))
      }
    }
  }
  return(model)
}
#------------------------------------------------------------------------------#
