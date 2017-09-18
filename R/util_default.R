#------------------------------------------------------------------------------#
#-------------- Some useful functions to fill in missing arguments  -----------#
#------------------------------------------------------------------------------#
#'@rdname set_default_algorithmic_parameters
#'@title set_default_algorithmic_parameters
#'@description Set the missing algorithmic parameters to their default values
#'@export
set_default_algorithmic_parameters = function(observations, model, algorithmic_parameters){
  # Number of theta-particles // Number of x-particles
  if (is.null(algorithmic_parameters$Ntheta)) {algorithmic_parameters$Ntheta = model$dimtheta*(2^7)}
  if (is.null(algorithmic_parameters$Nx)) {algorithmic_parameters$Nx = 2^ceiling(log2(ncol(observations)*model$dimY))}
  # Adaptive Nx or not // Maximum value for Nx // Acceptance rate threshold to trigger the increase Nx step
  if (is.null(algorithmic_parameters$adaptNx)) {algorithmic_parameters$adaptNx = TRUE}
  if (is.null(algorithmic_parameters$Nx_max)) {algorithmic_parameters$Nx_max = Inf}
  if (is.null(algorithmic_parameters$min_acceptance_rate)) {algorithmic_parameters$min_acceptance_rate = 0.20}
  # ESS threshold to trigger the resample-move steps // Number of moves per rejuvenation steps
  if (is.null(algorithmic_parameters$ess_threshold)) {algorithmic_parameters$ess_threshold = 0.5}
  if (is.null(algorithmic_parameters$nmoves)) {algorithmic_parameters$nmoves = 1}
  # Compute the Hyvarinen score or not
  if (is.null(algorithmic_parameters$hscore)) {algorithmic_parameters$hscore = TRUE}
  # For discrete observations, specify which differencing scheme to use
  # ("forward" or "central", default is set to "central")
  if (algorithmic_parameters$hscore && model$observation_type == "discrete") {
    if (is.null(algorithmic_parameters$discrete_diff_type)) {algorithmic_parameters$discrete_diff_type = "central"}
  }
  # Keep all the history of theta-particles // x-particles // byproducts (e.g. auxiliary Kalman filters)
  if (is.null(algorithmic_parameters$store_thetas_history)) {algorithmic_parameters$store_thetas_history = FALSE}
  if (is.null(algorithmic_parameters$store_X_history)) {algorithmic_parameters$store_X_history = FALSE}
  if (is.null(algorithmic_parameters$store_byproducts_history)) {algorithmic_parameters$store_byproducts_history = FALSE}
  # Keep the most recent theta-particles // x-particles // byproducts (e.g. auxiliary Kalman filters)
  if (is.null(algorithmic_parameters$store_last_thetas)) {algorithmic_parameters$store_last_thetas = TRUE}
  if (is.null(algorithmic_parameters$store_last_X)) {algorithmic_parameters$store_last_X = TRUE}
  if (is.null(algorithmic_parameters$store_last_byproducts)) {algorithmic_parameters$store_last_byproducts = TRUE}
  # Display progress bar // Display diagnostic at each step
  if (is.null(algorithmic_parameters$progress)) {algorithmic_parameters$progress = FALSE}
  if (is.null(algorithmic_parameters$verbose)) {algorithmic_parameters$verbose = FALSE}
  # Save intermediary results in RDS file // Allocate a time budget for the computation
  # WARNING: allocating a time budget is only relevant when the option "save" is on, otherwise
  # all the results will be lost if the time limit is reached and the computation gets interrupted.
  if (is.null(algorithmic_parameters$save)) {algorithmic_parameters$save = FALSE}
  # Save results once every save_stepsize observations (only relevant if save = TRUE)
  if (is.null(algorithmic_parameters$save_stepsize)) {algorithmic_parameters$save_stepsize = 1}
  if (algorithmic_parameters$save) {
    # Default schedule of times at which to save results
    if (is.null(algorithmic_parameters$save_schedule)) {
      algorithmic_parameters$save_schedule = seq(1,ncol(observations),algorithmic_parameters$save_stepsize)
      algorithmic_parameters$save_schedule = c(algorithmic_parameters$save_schedule, ncol(observations))
    }
  }
  if (is.null(algorithmic_parameters$time_budget)) {algorithmic_parameters$time_budget = NULL}
  # File name where the intermediary results will be saved
  # WARNING: must be an RDS file (.rds). Save in the working directory by default with timestamp as name.
  if (is.null(algorithmic_parameters$savefilename)) {
    algorithmic_parameters$savefilename = paste("results_",format(Sys.time(), "%Y-%m-%d_%H-%M-%S"),".rds",sep="")
  }
  # Resampling scheme: the default option is systematic resampling
  if (is.null(algorithmic_parameters$resampling)) {
    algorithmic_parameters$resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1))
  }
  # Proposal for rejuvenation steps: the default is independent draws from a fitted mixture of Normals
  if (is.null(algorithmic_parameters$proposalmove)) {
    algorithmic_parameters$proposalmove = get_proposal_mixture()
  }
  # Use kernel density estimators to compute the log-predictives and their derivatives
  if (is.null(algorithmic_parameters$use_dde)) {algorithmic_parameters$use_dde = FALSE}
  # parameters for density estimators via local regression
  if (algorithmic_parameters$use_dde) {
    if (is.null(algorithmic_parameters$dde_options)) {
      # options for density estimation
      algorithmic_parameters$dde_options = list(Ny = 10^4,
                                            sigma2_order0 = 0.001,
                                            sigma2_order1 = 0.01,
                                            sigma2_order2 = 0.01,
                                            nb_steps = Inf)
    } else {
      if (is.null(algorithmic_parameters$dde_options$Ny)) {algorithmic_parameters$dde_options$Ny = 10^4}
      if (is.null(algorithmic_parameters$dde_options$sigma2_order0)) {algorithmic_parameters$dde_options$sigma2_order0 = 0.001}
      if (is.null(algorithmic_parameters$dde_options$sigma2_order1)) {algorithmic_parameters$dde_options$sigma2_order1 = 0.002}
      if (is.null(algorithmic_parameters$dde_options$sigma2_order2)) {algorithmic_parameters$dde_options$sigma2_order2 = 0.01}
      if (is.null(algorithmic_parameters$dde_options$nb_steps)) {algorithmic_parameters$dde_options$nb_steps = Inf}
    }
  }
  # Reduce variance: if TRUE, generate Nc additional (temporary) particles
  # to compute the Hyvarinen score
  if (is.null(algorithmic_parameters$reduce_variance)) {algorithmic_parameters$reduce_variance = FALSE}
  if (algorithmic_parameters$reduce_variance) {
    if (is.null(algorithmic_parameters$Nc)) {algorithmic_parameters$Nc = 2*algorithmic_parameters$Ntheta}
    if (is.null(algorithmic_parameters$Ncx)) {algorithmic_parameters$Ncx = 2*algorithmic_parameters$Nx}
  }
  # return complete set of parameters
  return(algorithmic_parameters)
}
#------------------------------------------------------------------------------#
#'@rdname set_default_model
#'@title set_default_model
#'@description Set the missing field of a model to their default values
#'@export
set_default_model = function(model){
  # Define the one-step predicitve density of the observation at time t given all the past from time 1 to (t-1)
  # This is only relevant when only the likelihood is specified and available
  if (!is.null(model$likelihood)){
    if (is.null(model$dpredictive)){
      # Some models require keeping track of some byproducts (e.g. Kalman filter recursion)
      if (is.null(model$initialize_byproducts)){
        model$dpredictive = function(observations,t,theta,log = TRUE){
          if (t==1){
            return(model$likelihood(observations,t,theta,log))
          } else {
            if (log){return(model$likelihood(observations,t,theta,log)-model$likelihood(observations,t-1,theta,log))}
            else {return(model$likelihood(observations,t,theta,log)/model$likelihood(observations,t-1,theta,log))}
          }
        }
      } else {
        model$dpredictive = function(observations,t,theta,byproduct,log = TRUE){
          if (t==1){
            return(model$likelihood(observations,t,theta,byproduct,log))
          } else {
            if (log){return(model$likelihood(observations,t,theta,byproduct,log)-model$likelihood(observations,t-1,theta,byproduct,log))}
            else {return(model$likelihood(observations,t,theta,byproduct,log)/model$likelihood(observations,t-1,theta,byproduct,log))}
          }
        }
      }
    }
  }
  # Define the likelihood of the observations from time 1 to t
  # This is only relevant when only the predictive is specified and available
  if (!is.null(model$dpredictive)){
    if (is.null(model$likelihood)){
      # Some models require keeping track of some byproducts (e.g. Kalman filter recursion)
      if (is.null(model$initialize_byproducts)){
        model$likelihood = function(observations,t,theta,log = TRUE){
          ll = 0
          for (i in 1:t) {ll = ll + model$dpredictive(observations,i,theta,log = TRUE)}
          if (log) {return(ll)}
          else {return(exp(ll))}
        }
      } else {
        model$likelihood = function(observations,t,theta,byproduct,log = TRUE){
          ll = 0
          for (i in 1:t){ll = ll + model$dpredictive(observations,i,theta,byproduct,log = TRUE)}
          if (log) {return(ll)}
          else {return(exp(ll))}
        }
      }
    }
  }
  # If the observations are discrete, lower and upper bounds defining the support for each component
  # of the observations are required to compute the Hyvarinen score
  if (tolower(model$observation_type) == "discrete") {
    if (is.null(model$lower)){
      model$lower = rep(-Inf, model$dimY)
    }
    if (is.null(model$upper)){
      model$upper = rep(Inf, model$dimY)
    }
  }
  if (tolower(model$observation_type) == "continuous") {
    # Define the derivatives of the observation log-density via numerical derivation (cf. numDeriv)
    # The function is vectorized with respect to the states Xt (dimX by Nx), so that it outputs:
    # >> the jacobian (Nx by dimY matrix: each row is the transpose of the corresponding gradients row-wise)
    # >> the Hessian diagonals (Nx by dimY matrix: each row is the diagonal coeffs of the corresponding Hessian)
    if (is.null(model$derivativelogdobs)){
      model$derivativelogdobs = function(Yt,Xt,t,theta){
        logdobs = function(y){model$dobs(y,Xt,t,theta,log = TRUE)}
        derivatives = get_derivatives_from_genD(genD(logdobs,Yt)$D,model$dimY)
        return (list(jacobian = derivatives$jacobian, hessiandiag = derivatives$hessiandiag))
      }
    }
    # Define the derivatives of the predictive log-density via numerical derivation (cf. numDeriv)
    # The function outputs:
    # >> the transpose of the gradient (1 by dimY)
    # >> the Hessian diagonal coefficients (1 by dimY)
    if (is.null(model$derivativelogdpredictive)){
      if (is.null(model$initialize_byproducts)){
        model$derivativelogdpredictive = function(observations,t,theta){
          if (t==1){
            logpred = function(y) {model$dpredictive(y,t,theta,log = TRUE)}
          } else {
            logpred = function(y) {model$dpredictive(cbind(observations[,1:(t-1),drop=FALSE],y),t,theta,log = TRUE)}
          }
          derivatives = get_derivatives_from_genD(genD(logpred,observations[,t,drop=FALSE])$D,model$dimY)
          return (list(jacobian = derivatives$jacobian, hessiandiag = derivatives$hessiandiag))
        }
      } else {
        model$derivativelogdpredictive = function(observations,t,theta,byproduct){
          if (t==1){
            logpred = function(y) {model$dpredictive(y,t,theta,byproduct,log = TRUE)}
          } else {
            logpred = function(y) {model$dpredictive(cbind(observations[,1:(t-1),drop=FALSE],y),t,theta,byproduct,log = TRUE)}
          }
          derivatives = get_derivatives_from_genD(genD(logpred,observations[,t,drop=FALSE])$D,model$dimY)
          return (list(jacobian = derivatives$jacobian, hessiandiag = derivatives$hessiandiag))
        }
      }
    }
  }
  return(model)
}
#------------------------------------------------------------------------------#

