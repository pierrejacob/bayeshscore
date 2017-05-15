#'@export bootstrap_particle_filter
#'@title bootstrap_particle_filter
#'@description This function implements the bootstrap particle filter. It returns log-evidence estimates, particles, and normalized weights.
#'@export
bootstrap_particle_filter <- function(observations, model, theta, algorithmic_parameters){
  Nx <- algorithmic_parameters$Nx
  resampling <- algorithmic_parameters$resampling
  nobservations <- ncol(observations)
  log_p_y_hat <- 0 #initialize estimate of p_theta(y_{1:nobservations}))
  X_history = list() # successive sets of Nx particles column-wise
  weight_history = list() # successive sets of normalized weights
  X = matrix(NA,nrow = model$dimX, ncol = Nx) # matrix of Nx particles column-wise (most recent)
  normW = rep(1/Nx, Nx) #vector of normalized weights
  X <- model$rinitial(theta,Nx) #initial step 1
  logW <- model$dobs(observations[,1],X,1,theta)
  maxlogW <- max(logW)
  W <- exp(logW - maxlogW)
  log_p_y_hat <- log_p_y_hat + log(mean(W)) + maxlogW #udpate likelihood estimate
  normW <- W / sum(W)
  X_history[[1]] = X
  weight_history[[1]] = normW
  #iterate for n = 2, ... T
  if (nobservations > 1){
    for (t in 2:nobservations) {
      ancestors <- resampling(normW) #sample the ancestors' indexes
      X <- X[,ancestors]
      if (is.null(dim(X))) X <- matrix(X, nrow = model$dimX)
      X <- model$rtransition(X, t, theta)
      logW <- model$dobs(observations[,t], X, t, theta)
      maxlogW <- max(logW)
      W <- exp(logW - maxlogW)
      log_p_y_hat <- log_p_y_hat + log(mean(W)) + maxlogW #udpate likelihood estimate
      normW <- W / sum(W)
      X_history[[t]] = X
      weight_history[[t]] = normW
    }
  }
  return(list(log_p_y_hat = log_p_y_hat, X_history = X_history, weight_history = weight_history))
}
