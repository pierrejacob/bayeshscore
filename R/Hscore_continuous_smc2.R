# This function computes successive prequential Hyvarinen score for continuous observations by running SMC^2.
# It also computes the successive log-evidence as a by-product.
hscore_continuous_smc2 <- function(observations, model, algorithmic_parameters){
  # Set default values for the missing fields
  algorithmic_parameters = set_default_algorithmic_parameters(algorithmic_parameters)
  model = set_default_model(model)
  # Parse algorithmic parameters and set flags accordingly
  nobservations = ncol(observations)
  Ntheta = algorithmic_parameters$Ntheta
  Nx = algorithmic_parameters$Nx
  nmoves = algorithmic_parameters$nmoves
  resampling = algorithmic_parameters$resampling
  adaptNx = algorithmic_parameters$adaptNx
  min_acceptance_rate = algorithmic_parameters$min_acceptance_rate
  ess_objective = algorithmic_parameters$ess_threshold*algorithmic_parameters$Ntheta
  # Monitor progress if needed
  if (algorithmic_parameters$progress) {
    print(paste("Started at:",Sys.time()))
    time_start = proc.time()
    progbar = txtProgressBar(min = 0,max = nobservations,style=3)
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

  ########## if we start from a proposal instead of the prior (e.g. improper prior)
  ########## then the weights should be initialized differently:
  ########## log(prior_density) - log(proposal_density) ???

  PFs = list() # list of particle filters (one for each theta)
  thetas_history[[1]] = thetas
  normw_history[[1]] = normw
  # Initialize filters (first observation passed as argument just to initialize the fields of PF)
  for (itheta in 1:Ntheta){
    theta = thetas[,itheta]
    PFs[[itheta]] = conditional_particle_filter(matrix(observations[,1],ncol = 1), model, theta, Nx)
    #the CPF performs a regular PF when no conditioning path is provided
  }
  # Assimilate observations one by one
  for (t in 1:nobservations){
    results = assimilate_one_smc2(thetas, PFs, t, observations, model, Ntheta, ess_objective,
                                  nmoves, resampling, logtargetdensities, logw, normw,
                                  algorithmic_parameters$verbose, adaptNx, min_acceptance_rate)
    thetas = results$thetas
    normw = results$normw
    logw = results$logw
    PFs = results$PFs
    logtargetdensities = results$logtargetdensities
    logevidence[t] = results$logcst
    # compute prequential H score here
    Hscore[t] = hincrementContinuous_smc2(t, model, observations[,t,drop=FALSE], thetas, normw, PFs, Ntheta)
    # do some book-keeping
    thetas_history[[t+1]] = thetas
    normw_history[[t+1]] = normw
    rejuvenation_times[t] = results$rejuvenation_time #successive times where resampling is triggered
    rejuvenation_accept_rate[t] = results$rejuvenation_accept_rate #successive acceptance rates
    increase_Nx_times[t] = results$increase_Nx_times #successive times where adaptation regarding Nx is triggered
    increase_Nx_values[t] = results$increase_Nx_values #successive values of Nx
    # Update progress bar if needed
    if (algorithmic_parameters$progress) {
      setTxtProgressBar(progbar, t)
    }
  }
  # Update progress bar if needed
  if (algorithmic_parameters$progress) {
    close(progbar)
    time_end = proc.time()-time_start
    cat(paste("Hscore: T = ",toString(nobservations),", Ntheta = ",toString(Ntheta),", Nx = ",toString(Nx),"\n",sep = ""))
    print(time_end)
  }
  return(list(thetas_history = thetas_history, normw_history = normw_history, logevidence = cumsum(logevidence),
              logtargetdensities = logtargetdensities, hscore = cumsum(Hscore),
              rejuvenation_times = rejuvenation_times[!is.na(rejuvenation_times)],
              rejuvenation_accept_rate = rejuvenation_accept_rate[!is.na(rejuvenation_accept_rate)],
              increase_Nx_times = increase_Nx_times[!is.na(increase_Nx_times)],
              increase_Nx_values = increase_Nx_values[!is.na(increase_Nx_values)]))
}
