#'@rdname smc
#'@title smc
#'@description This function runs the SMC algorithm, using tempering.
#'It also computes the log-evidence and the prequential Hyvarinen score (optional).
#'@export
smc = function(observations, model, algorithmic_parameters){
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
  incr_logevidence = array(NA,dim = c(nobservations)) #incremental log-evidence at successive times t
  if (algorithmic_parameters$hscore) {incr_hscore = array(NA,dim = c(nobservations))} # OPTIONAL: incremental Hyvarinen score at successive times t
  rejuvenation_times = c() #successive times where resampling is triggered
  rejuvenation_rate = c() #successive acceptance rates of resampling
  thetas_history = list() #successive sets of particles theta
  normw_history = list() #successive sets of normalized weights for theta
  # if (is.null(algorithmic_parameters$rinitial_theta)){
  thetas = model$rprior(Ntheta)
  # } else {
  # thetas = algorithmic_parameters$rinitial_theta(Ntheta)
  # }
  logtargetdensities = apply(thetas, 2, model$dprior) # log target density evaluations at current particles
  normw = rep(1/Ntheta, Ntheta) # normalized weights
  logw = rep(0, Ntheta) # log normalized weights

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
    results = assimilate_one_smc(thetas,byproducts,t,observations,model,logtargetdensities,logw,normw,algorithmic_parameters)
    # Update the particles theta and compute the log-evidence
    thetas = results$thetas
    normw = results$normw
    logw = results$logw
    logtargetdensities = results$logtargetdensities
    incr_logevidence[t] = results$logcst
    if (!is.null(byproducts)) {byproducts = results$byproducts}
    # OPTIONAL: compute the incremental hscore for continuous observations
    if (algorithmic_parameters$hscore && (observation_type=="continuous")) {
      incr_hscore[t] = hincrementContinuous_smc(t, model, observations,thetas,normw,byproducts,Ntheta)
    }
    # do some book-keeping
    ESS[t] = results$ESS
    if (!is.na(results$rejuvenation_time)) {rejuvenation_times = c(rejuvenation_times, results$rejuvenation_time)}
    if (!is.na(results$rejuvenation_rate)) {rejuvenation_rate = c(rejuvenation_rate, results$rejuvenation_rate)}
    if (algorithmic_parameters$store_theta){
      thetas_history[[t+1]] = thetas
      normw_history[[t+1]] = normw
    }
    # Update progress bar if needed
    if (algorithmic_parameters$progress) {
      setTxtProgressBar(progbar, t)
    }
    # save partial results if needed
    if (algorithmic_parameters$save) {
      # save the variables required to resume and proceed further, in case of interrupted run
      # NOTE: thetas and normw, are retrievable from the history stored in results_so_far
      required_to_resume = list(t = t, logw = logw, logtargetdensities = logtargetdensities,
                                byproducts = byproducts, observations = observations, model = model,
                                algorithmic_parameters = algorithmic_parameters)
      # save the results obtained up to this time
      results_so_far = list(thetas_history = thetas_history, normw_history = normw_history,
                            incr_logevidence = incr_logevidence[1:t], incr_hscore = incr_hscore[1:t],  ESS = ESS[1:t],
                            rejuvenation_times = rejuvenation_times, rejuvenation_rate = rejuvenation_rate)
      # if the history of theta-particles is not saved, just keep the most recent ones
      if (!algorithmic_parameters$store_theta){
        required_to_resume$thetas = thetas
        required_to_resume$normw = normw
      }
      # save into RDS file
      saveRDS(c(required_to_resume,results_so_far),file = algorithmic_parameters$savefilename)
    }
  }
  # Update progress bar if needed
  if (algorithmic_parameters$progress) {
    close(progbar)
    time_end = proc.time()-time_start
    cat(paste("SMC: T = ",toString(nobservations),", Ntheta = ",toString(Ntheta),"\n",sep=""))
    print(time_end)
  }
  return (list(thetas_history = thetas_history, normw_history = normw_history, logevidence = cumsum(incr_logevidence),
               logtargetdensities = logtargetdensities, hscore = cumsum(incr_hscore), ESS = ESS,
               rejuvenation_times = rejuvenation_times, rejuvenation_rate = rejuvenation_rate))
}
