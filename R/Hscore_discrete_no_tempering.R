#'@rdname hscore_discrete_no_tempering
#'@title hscore_discrete_no_tempering
#'@description This function computes successive prequential Hyvarinen score for discrete observations by running the smc^2 algorithm. It also computes the successive log-evidence as a by-product.
#'@export
hscore_discrete_no_tempering <- function(observations, model, algorithmic_parameters){
  # Extract algorithmic parameters and set flags accordingly
  Ntheta <- algorithmic_parameters$Ntheta
  nobservations <- ncol(observations)
  # If no number of X-particles Nx is specified, use adaptive Nx starting with Nx = 128
  if (is.null(algorithmic_parameters$Nx)) {
    adaptNx = TRUE
    Nx = 2^7
    algorithmic_parameters$Nx = Nx
    # Trigger increase in Nx when acceptance rate drops below threshold (default is 10%)
    if (is.null(algorithmic_parameters$min_acceptance_rate)) {
      min_acceptance_rate = 0.10
    }
    else {
      min_acceptance_rate = algorithmic_parameters$min_acceptance_rate
    }
  }
  else {
    adaptNx = FALSE
    Nx = algorithmic_parameters$Nx
  }
  if (is.null(algorithmic_parameters$progress)) algorithmic_parameters$progress <- FALSE
  if (is.null(algorithmic_parameters$store)) algorithmic_parameters$store <- FALSE
  if (algorithmic_parameters$progress) {
    print(paste("Started at:",Sys.time()))
    progbar = txtProgressBar(min = 0,max = nobservations,style=3)
    count = 0
    time_start = proc.time()
  }
  # Initialize empty arrays and lists to store the results
  ESS = array(NA,dim = c(nobservations)) #ESS at successive times t
  Hscore = array(NA,dim = c(nobservations)) #prequential Hyvarinen score at successive times t
  logevidence = array(NA,dim = c(nobservations)) #log-evidence at successive times t
  rejuvenation_times <- c() #successive times where resampling is triggered
  rejuvenation_accept_rate <- c() #successive acceptance rates of resampling
  increase_Nx_times <- c() #successive times where adaptation regarding Nx is triggered
  increase_Nx_values <- c() #successive values of Nx
  # Initialize array of particles X targeting predictive distributions (most recent)
  Xpred = array(NA,dim = c(Nx, model$dimX, Ntheta))
  # Initialize SMC2 by sampling from prior (default), or from specified proposal
  # (e.g. relevant in case the prior is vague)
  if (is.null(algorithmic_parameters$rinitial_theta)){
    thetas <- model$rprior(Ntheta)
  }
  else {
    thetas <- algorithmic_parameters$rinitial_theta(Ntheta)
  }
  # Construct list of trees to store paths (one tree for each theta)
  trees = list()
  for (itheta in 1:Ntheta){
    trees[[itheta]] = new(TreeClass, Nx, 10*Nx, model$dimX)
  }
  # Initialize particle filter + one step (i.e. including 1 observations)
  first_step <- filter_first_step(observations[,1], model, thetas, trees, algorithmic_parameters)
  X = first_step$X  #Nx particles (most recent) for each theta (size = Nx,dimX,Ntheta)
  xnormW = first_step$xnormW  #matrix of corresponding normalized X-weights (size = Nx,Ntheta)
  Xpred = X
  XprednormW = matrix(1/Nx, nrow = Nx, ncol = Ntheta) #matrix of normalized weights for Xpred
  log_z = first_step$log_z  #matrix of log-likelihood estimates (size = Ntheta)
  trees = first_step$trees #list of trees to extract X-paths
  # Initialize the weights for theta
  if (is.null(algorithmic_parameters$rinitial_theta)){
    # Equal initial weights when sampling from the prior initially (default)
    thetalogw <- rep(0,Ntheta)
    thetanormw <- rep(1/Ntheta,Ntheta)
  }
  else {
    # if we initialize the first theta by using a proposal different from the prior, the first initial
    # weights are not all equal but are importance weights targeting the prior. This happens before
    # seeing the first observation
    thetalogw <- rep(NA,Ntheta)
    for (m in 1:Ntheta) {
      thetalogw[m] <- model$dprior(thetas[m,], log = TRUE) - algorithmic_parameters$dinitial_theta(thetas[m,], log = TRUE)
    }
    maxlogW <- max(thetalogw) #avoids overflow when exponentiating
    W <- exp(thetalogw - maxlogW) #computes actual unnormalized weights for theta
    thetanormw <- W / sum(W) #normalize weights for theta
  }
  maxlogz = max(log_z) #avoids overflow when exponentiating
  actual_z = exp(log_z - maxlogz) #actual z up to a multiplicative constant
  # compute prequential H score and log-evidence
  Hscore[1] = SHd(1,model,observations[,1],thetas,thetanormw,Xpred,XprednormW,Ntheta)
  logevidence[1] = log(sum(actual_z*thetanormw)) + maxlogz
  # update the weights after seeing the first observation
  thetalogw <- thetalogw + log_z #update log-weights for theta
  maxlogW <- max(thetalogw) #avoids overflow when exponentiating
  W <- exp(thetalogw - maxlogW) #computes actual unnormalized weights for theta
  thetanormw <- W / sum(W) #normalize weights for theta
  # Store thetas and weights if needed
  thetas_history <- list()
  weights_history <- list()
  if (algorithmic_parameters$store){
    thetas_history[[1]] <- thetas
    weights_history[[1]] <- thetanormw
  }
  # Compute ESS and resample if needed
  ESS[1] <- getESS(thetanormw)
  if ((ESS[1]/Ntheta) < 0.5){
    rejuvenation = rejuvenation_step(observations, 1, model, thetas, thetanormw, X, xnormW, log_z, trees, algorithmic_parameters)
    thetas = rejuvenation$thetas
    thetalogw = rep(log(1/Ntheta),Ntheta)
    X = rejuvenation$X
    xnormW = rejuvenation$xnormW
    log_z = rejuvenation$log_z
    trees = rejuvenation$trees
    rejuvenation_times <- c(rejuvenation_times, 1)
    rejuvenation_accept_rate <- c(rejuvenation_accept_rate, rejuvenation$accept_rate)
    if (adaptNx){
      if (rejuvenation$accept_rate < min_acceptance_rate){
        # Increase the number Nx of particles for each theta
        new_filter = increase_Nx_no_tempering(observations, 1, model, thetas, xnormW, trees, algorithmic_parameters)
        X = new_filter$X
        xnormW = new_filter$xnormW
        log_z = new_filter$log_z
        trees = new_filter$trees
        Nx = new_filter$new_Nx
        algorithmic_parameters$Nx = Nx
        increase_Nx_times <- c(increase_Nx_times,1)
        increase_Nx_values <- c(increase_Nx_values,new_filter$new_Nx)
      }
    }
  }
  # Update progress bar if needed
  if (algorithmic_parameters$progress) {
    count = count + 1
    setTxtProgressBar(progbar, count)
  }
  # Iterate over the next observations
  for (t in 2:nobservations){
    Xpred = filter_predict(t, model, thetas, X, algorithmic_parameters)
    XprednormW = xnormW
    # compute prequential H score
    Hscore[t] = Hscore[t-1] + SHd(t,model,observations[,t],thetas,thetanormw,Xpred,XprednormW,Ntheta)
    # Propagate X-particles
    next_step <- filter_next_step(observations[,t], t, model, thetas, X, xnormW, trees, algorithmic_parameters)
    X = next_step$X
    xnormW = next_step$xnormW
    log_z_incremental = next_step$log_z_incremental
    trees = next_step$trees
    log_z <- log_z + log_z_incremental
    # update log-evidence estimator
    maxlogz_incremental = max(log_z_incremental)
    logevidence[t] <- logevidence[t-1] + log(sum(exp(log_z_incremental - maxlogz_incremental)*thetanormw)) + maxlogz_incremental
    # Reweighting
    thetalogw <- thetalogw + log_z_incremental #update log-weights for theta
    maxlogW <- max(thetalogw) #avoids overflow when exponentiating
    W <- exp(thetalogw - maxlogW) #computes actual unnormalized weights for theta
    thetanormw <- W / sum(W) #normalize weights for theta
    # Store thetas and weights if needed
    if (algorithmic_parameters$store){
      thetas_history[[t]] <- thetas
      weights_history[[t]] <- thetanormw
    }
    # Compute ESS and resample if needed
    ESS[t] <- getESS(thetanormw)
    if ((ESS[t]/Ntheta) < 0.5){
      rejuvenation = rejuvenation_step(observations, t, model, thetas, thetanormw, X, xnormW, log_z, trees, algorithmic_parameters)
      thetas = rejuvenation$thetas
      thetalogw = rep(log(1/Ntheta),Ntheta)
      X = rejuvenation$X
      xnormW = rejuvenation$xnormW
      log_z = rejuvenation$log_z
      trees = rejuvenation$trees
      rejuvenation_times <- c(rejuvenation_times, t)
      rejuvenation_accept_rate <- c(rejuvenation_accept_rate, rejuvenation$accept_rate)
      if (adaptNx){
        if (rejuvenation$accept_rate < min_acceptance_rate){
          # Increase the number Nx of particles for each theta
          new_filter = increase_Nx_no_tempering(observations, t, model, thetas, xnormW, trees, algorithmic_parameters)
          X = new_filter$X
          xnormW = new_filter$xnormW
          log_z = new_filter$log_z
          trees = new_filter$trees
          Nx = new_filter$new_Nx
          algorithmic_parameters$Nx = Nx
          increase_Nx_times <- c(increase_Nx_times,t)
          increase_Nx_values <- c(increase_Nx_values,new_filter$new_Nx)
        }
      }
    }
    # Update progress bar if needed
    if (algorithmic_parameters$progress) {
      count = count + 1
      setTxtProgressBar(progbar, count)
    }
  }
  # Update progress bar if needed
  if (algorithmic_parameters$progress) {
    close(progbar)
    time_end = proc.time()-time_start
    cat(paste("Hscore: T = ",toString(nobservations),", Ntheta = ",toString(Ntheta),
              ", Nx = ",toString(Nx),"\n",sep = ""))
    print(time_end)
  }
  # Return results
  if (algorithmic_parameters$store){
    return (list(hscore = Hscore, logevidence = logevidence, ESS = ESS, thetas = thetas, thetanormw = thetanormw,
                 thetas_history = thetas_history, weights_history = weights_history, trees = trees, xnormW = xnormW,
                 rejuvenation_times = rejuvenation_times, rejuvenation_accept_rate = rejuvenation_accept_rate,
                 increase_Nx_times = increase_Nx_times, increase_Nx_values = increase_Nx_values))
  }
  else {
    return (list(hscore = Hscore, logevidence = logevidence, ESS = ESS, thetas = thetas, thetanormw = thetanormw,
                 trees = trees, xnormW = xnormW, rejuvenation_times = rejuvenation_times,
                 rejuvenation_accept_rate = rejuvenation_accept_rate,
                 increase_Nx_times = increase_Nx_times, increase_Nx_values = increase_Nx_values))
  }
}

