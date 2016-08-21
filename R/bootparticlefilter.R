#'@export

bootstrap_particle_filter <- function(observations, model, theta, algorithmic_parameters){
  Nx <- algorithmic_parameters$Nx
  resampling <- algorithmic_parameters$resampling
  nobservations <- ncol(observations)
  log_p_y_hat <- 0 #initialize estimate of p_theta(y_{1:nobservations}))
  X = matrix(NA,nrow = Nx, ncol = model$dimX) # matrix of Nx particles row-wise
  normW = rep(1/Nx, Nx) #vector of normalized weights
  X <- model$rinitial(theta,Nx) #initial step 1
  logW <- model$dobs(observations[,1],X,1,theta)
  maxlogW <- max(logW)
  W <- exp(logW - maxlogW)
  log_p_y_hat <- log_p_y_hat + log(mean(W)) + maxlogW #udpate likelihood estimate
  normW <- W / sum(W)
  #iterate for n = 2, ... T
  for (t in 2:nobservations) {
    ancestors <- resampling(normW) #sample the ancestors' indexes
    X <- X[ancestors,]
    if (is.null(dim(X))) X <- matrix(X, ncol = model$dimX)
    X <- model$rtransition(X, t, theta)
    logW <- model$dobs(observations[,t], X, t, theta)
    maxlogW <- max(logW)
    W <- exp(logW - maxlogW)
    log_p_y_hat <- log_p_y_hat + log(mean(W)) + maxlogW #udpate likelihood estimate
    normW <- W / sum(W)
  }
  return(list(log_p_y_hat = log_p_y_hat, X = X, normW = normW))
}
