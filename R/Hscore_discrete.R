#'@rdname hscore_discrete
#'@title hscore_discrete
#'@description This function computes successive prequential Hyvarinen score for discrete observations by running the smc^2 algorithm. It also computes the successive log-evidence as a by-product.
#'@export
hscore_discrete <- function(observations, model, algorithmic_parameters){
  # Extract algorithmic parameters and set flags accordingly
  Ntheta = algorithmic_parameters$Ntheta
  nmoves = algorithmic_parameters$nmoves
  resampling = algorithmic_parameters$resampling
  ess_objective = algorithmic_parameters$ess_threshold*algorithmic_parameters$Ntheta
  nobservations = ncol(observations)
  # If no number of X-particles Nx is specified, use adaptive Nx starting with Nx = 128
  if (is.null(algorithmic_parameters$Nx)) {
    adaptNx = TRUE
    Nx = 2^7
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
  # Set flags for progress bar and storage of particle history
  if (is.null(algorithmic_parameters$progress)) algorithmic_parameters$progress = FALSE
  if (is.null(algorithmic_parameters$store)) algorithmic_parameters$store = FALSE
  if (algorithmic_parameters$progress) {
    print(paste("Started at:",Sys.time()))
    time_start = proc.time()
  }
  # Initialize empty arrays and lists to store the results
  ESS = array(NA,dim = c(nobservations)) #ESS at successive times t
  Hscore = array(NA,dim = c(nobservations)) #prequential Hyvarinen score at successive times t
  logevidence = array(NA,dim = c(nobservations)) #log-evidence at successive times t
  rejuvenation_times = array(NA,dim = c(nobservations)) #successive times where resampling is triggered
  rejuvenation_accept_rate = array(NA,dim = c(nobservations)) #successive acceptance rates of resampling
  increase_Nx_times = array(NA,dim = c(nobservations)) #successive times where adaptation regarding Nx is triggered
  increase_Nx_values = array(NA,dim = c(nobservations)) #successive values of Nx
  thetas_history = list() #successive sets of particles theta
  normw_history = list() #successive sets of normalized weights for theta
  #
  # # Initialize SMC2 by sampling from prior (default), or from specified proposal
  # # (e.g. relevant in case the prior is vague)
  # if (is.null(algorithmic_parameters$rinitial_theta)){
  thetas = model$rprior(Ntheta)
  # } else {
  # thetas = algorithmic_parameters$rinitial_theta(Ntheta)
  # }
  #
  logtargetdensities = apply(thetas, 1, model$dprior) # log target density evaluations at current particles
  normw = rep(1/Ntheta, Ntheta) # normalized weights
  logw = rep(0, Ntheta) # log normalized weights
  PFs = list() # list of particle filters (one for each theta)
  thetas_history[[1]] = thetas
  normw_history[[1]] = normw
  # Initialize array of particles X targeting predictive distributions (most recent)
  Xprevious = array(NA,dim = c(Nx, model$dimX, Ntheta))
  Xpred = array(NA,dim = c(Nx, model$dimX, Ntheta))
  XnormW_previous = matrix(1/Nx, nrow = Nx, ncol = Ntheta) #matrix of normalized weights for X at previous step
  # Initialize filters (first observation passed as argument just to initialize the fields of PF)
  for (itheta in 1:Ntheta){
    theta = thetas[itheta,]
    PFs[[itheta]] = conditional_particle_filter(matrix(observations[,1],ncol = 1), model, theta, Nx)
    #the CPF performs a regular PF when no conditioning path is provided
    Xpred[,,itheta] = PFs[[itheta]]$X
  }
  # Assimilate observations one by one
  for (t in 1:nobservations){
    # Construct particles targeting the one-step-ahead predictive (need to reconstruct since size Nx might change)
    if (t > 1){
      Nx = PFs[[1]]$Nx
      Xpred = array(NA,dim = c(Nx, model$dimX, Ntheta)) # (need to reconstruct since size Nx might change)
      for (itheta in 1:Ntheta){
        X = Xprevious[,,itheta]
        if (is.null(dim(X))){
          Xpred[,,itheta] = model$rtransition(matrix(X,nrow = Nx), t, thetas[itheta,])
        }
        else{
          Xpred[,,itheta] = model$rtransition(X, t, thetas[itheta,])
        }
      }
    }
    # compute prequential H score (with theta from time t-1, see formula in the paper)
    Hscore[t] = SHd(t,model,observations[,t],thetas,normw,Xpred,XnormW_previous,Ntheta)
    # assimilate the next observation
    results = assimilate_next(thetas, PFs, t, observations, model, Ntheta, ess_objective,
                              nmoves, resampling, logtargetdensities, logw, normw,
                              algorithmic_parameters$progress, adaptNx, min_acceptance_rate)
    # Update the particles theta and compute the log-evidence
    thetas = results$thetas
    normw = results$normw
    logw = results$logw
    PFs = results$PFs
    logtargetdensities = results$logtargetdensities
    logevidence[t] = results$logcst
    #matrix of normalized weights for X at previous step (need to reconstruct since size Nx might change)
    Nx = PFs[[1]]$Nx
    Xprevious = array(NA,dim = c(Nx, model$dimX, Ntheta))
    XnormW_previous = matrix(NA, nrow = PFs[[1]]$Nx, ncol = Ntheta)
    for (itheta in 1:Ntheta){
      Xprevious[,,itheta] = PFs[[itheta]]$X
      XnormW_previous[,itheta] = PFs[[itheta]]$xnormW
    }
    # do some book-keeping
    thetas_history[[t+1]] = thetas
    normw_history[[t+1]] = normw
    rejuvenation_times[t] = results$rejuvenation_time #successive times where resampling is triggered
    rejuvenation_accept_rate[t] = results$rejuvenation_accept_rate #successive acceptance rates
    increase_Nx_times[t] = results$increase_Nx_times #successive times where adaptation regarding Nx is triggered
    increase_Nx_values[t] = results$increase_Nx_values #successive values of Nx
  }
  return(list(thetas_history = thetas_history, normw_history = normw_history, logevidence = cumsum(logevidence),
              logtargetdensities = logtargetdensities, Hscore = cumsum(Hscore),
              rejuvenation_times = rejuvenation_times[!is.na(rejuvenation_times)],
              rejuvenation_accept_rate = rejuvenation_accept_rate[!is.na(rejuvenation_accept_rate)],
              increase_Nx_times = increase_Nx_times[!is.na(increase_Nx_times)],
              increase_Nx_values = increase_Nx_values[!is.na(increase_Nx_values)]))
}
