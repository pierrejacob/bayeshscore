#'@rdname conditional_particle_filter
#'@title conditional_particle_filter
#'@description This function implements the conditional particle filter. If no conditioning path is provided, it performs a regular bootstrap particle filter. It returns log-evidence estimates, particles, normalized weights, tree of paths, and one sample final path.
#'@export
conditional_particle_filter <- function(observations, model, theta, algorithmic_parameters, path = NULL){
  # TreeClass = algorithmic_parameters$TreeClass
  Tree <- new(TreeClass, algorithmic_parameters$Nx, 10*algorithmic_parameters$Nx, model$dimX)
  # we expect a path as a matrix of size dimX x nobservations
  cpf <- TRUE
  if (is.null(path)){
    cpf <- FALSE
  }
  Nx <- algorithmic_parameters$Nx
  nobservations <- ncol(observations)
  log_p_y_hat <- 0 #initialize estimate of log-evidence
  X = matrix(NA,nrow = Nx, ncol = model$dimX) # matrix of Nx particles row-wise
  normW = rep(1/Nx, Nx) #vector of normalized weights
  X <- model$rinitial(theta,Nx) #initial step 1
  logW <- model$dobs(observations[,1],X,1,theta)
  maxlogW <- max(logW)
  W <- exp(logW - maxlogW)
  log_p_y_hat <- log_p_y_hat + log(mean(W)) + maxlogW #udpate likelihood estimate
  normW <- W / sum(W)
  #
  if (cpf){
    X[Nx,] <- path[,1]
  }
  Tree$init(t(X))
  #iterate for n = 2, ... T
  if (nobservations > 1){
    for (t in 2:nobservations) {
      ancestors <- sample(x = 1:Nx, size = Nx, replace = TRUE, prob = normW) #sample the ancestors' indexes
      if (cpf){
        ancestors[Nx] <- Nx
      }
      X <- X[ancestors,]
      if (is.null(dim(X))) X <- matrix(X, ncol = model$dimX)
      X <- model$rtransition(X, t, theta)
      if (cpf){
        X[Nx,] <- path[,t]
      }
      logW <- model$dobs(observations[,t], X, t, theta)
      maxlogW <- max(logW)
      W <- exp(logW - maxlogW)
      log_p_y_hat <- log_p_y_hat + log(mean(W)) + maxlogW #udpate likelihood estimate
      normW <- W / sum(W)
      Tree$update(t(X), ancestors-1)
    }
  }
  new_path <- Tree$get_path(sample(x = 0:(Nx-1), size = 1, replace = TRUE, prob = normW))
  return(list(log_p_y_hat = log_p_y_hat, X = X, normW = normW, tree = Tree, path = new_path))
}
