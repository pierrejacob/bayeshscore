#'@rdname smc
#'@title smc
#'@description This function runs the SMC algorithm, using tempering.
#'It also computes the log-evidence and the prequential Hyvarinen score (optional).
#'@export
smc <- function(observations, model, algorithmic_parameters){
  # Parse algorithmic parameters and set flags accordingly
  nobservations = ncol(observations)
  Ntheta = algorithmic_parameters$Ntheta
  nmoves = algorithmic_parameters$nmoves
  resampling = algorithmic_parameters$resampling
  ess_objective = algorithmic_parameters$ess_threshold*algorithmic_parameters$Ntheta
  if (algorithmic_parameters$hscore) {observation_type = tolower(model$observation_type)}
  # Monitor progress if needed
  if (algorithmic_parameters$progress) {
    print(paste("Started at:",Sys.time()))
    progbar = txtProgressBar(min = 0,max = nobservations,style=3)
    time_start = proc.time()
  }
  # Initialize empty arrays and lists to store the results
  ESS = array(NA,dim = c(nobservations)) #ESS at successive times t
  logevidence = array(NA,dim = c(nobservations)) #log-evidence at successive times t
  if (algorithmic_parameters$hscore) {incr_hscore = array(NA,dim = c(nobservations))} # OPTIONAL: incremental Hyvarinen score at successive times t
  rejuvenation_times = array(NA,dim = c(nobservations)) #successive times where resampling is triggered
  rejuvenation_accept_rate = array(NA,dim = c(nobservations)) #successive acceptance rates of resampling
  thetas_history = list() #successive sets of particles theta
  normw_history = list() #successive sets of normalized weights for theta
  # if (is.null(algorithmic_parameters$rinitial_theta)){
  thetas <- model$rprior(Ntheta)
  # } else {
  # thetas <- algorithmic_parameters$rinitial_theta(Ntheta)
  # }
  logtargetdensities <- apply(thetas, 2, model$dprior) # log target density evaluations at current particles
  normw <- rep(1/Ntheta, Ntheta) # normalized weights
  logw <- rep(0, Ntheta) # log normalized weights

  ########## if we start from a proposal instead of the prior (e.g. improper prior)
  ########## then the weights should be initialized differently:
  ########## log(prior_density) - log(proposal_density) ???
  if (algorithmic_parameters$store_theta){
    thetas_history[[1]] = thetas
    normw_history[[1]] = normw
  }
  # initialize possible byproducts (e.g. Kalman filters, etc ...)
  byproducts = list()
  if (!is.null(model$initialize_byproducts)) {
    for (itheta in 1:Ntheta){
      byproducts[[itheta]] = model$initialize_byproducts(thetas[,itheta], observations)
    }
  } else {
    byproducts = NULL
  }
  # Assimilate observations one by one
  for (t in 1:nobservations){
    # OPTIONAL: compute the incremental hscore for discrete observations
    if (algorithmic_parameters$hscore && (observation_type=="discrete")) {
      # compute incremental H score (with theta from time t-1, see formula in the paper)
      incr_hscore[t] = Hd_smc(t,model,observations,thetas,normw,Ntheta,byproducts)
    }
    # Assimilate the next observation
    results <- assimilate_one_smc(thetas, t, observations, model, Ntheta, ess_objective,
                                  nmoves, resampling, logtargetdensities, logw, normw,
                                  byproducts, algorithmic_parameters$verbose)
    # Update the particles theta and compute the log-evidence
    thetas <- results$thetas
    normw <- results$normw
    logw <- results$logw
    logtargetdensities <- results$logtargetdensities
    logevidence[t] <- results$logcst
    if (!is.null(byproducts)) {byproducts = results$byproducts}
    # OPTIONAL: compute the incremental hscore for continuous observations
    if (algorithmic_parameters$hscore && (observation_type=="continuous")) {
      incr_hscore[t] = hincrementContinuous_smc(t, model, observations,thetas,normw,byproducts,Ntheta)
    }
    # do some book-keeping
    rejuvenation_times[t] = results$rejuvenation_time #successive times where resampling is triggered
    rejuvenation_accept_rate[t] = results$rejuvenation_accept_rate #successive acceptance rates
    if (algorithmic_parameters$store_theta){
      thetas_history[[t+1]] = thetas
      normw_history[[t+1]] = normw
    }
    # Update progress bar if needed
    if (algorithmic_parameters$progress) {
      setTxtProgressBar(progbar, t)
    }
  }
  # Update progress bar if needed
  if (algorithmic_parameters$progress) {
    close(progbar)
    time_end = proc.time()-time_start
    cat(paste("SMC: T = ",toString(nobservations),", Ntheta = ",toString(Ntheta),"\n",sep=""))
    print(time_end)
  }
  return (list(thetas_history = thetas_history, normw_history = normw_history, logevidence = cumsum(logevidence),
               logtargetdensities = logtargetdensities, hscore = cumsum(incr_hscore),
               rejuvenation_times = rejuvenation_times[!is.na(rejuvenation_times)],
               rejuvenation_accept_rate = rejuvenation_accept_rate[!is.na(rejuvenation_accept_rate)]))
}

