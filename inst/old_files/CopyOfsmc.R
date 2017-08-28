#'@rdname smc
#'@title smc
#'@description This function runs the SMC algorithm, using tempering.
#'It also computes the log-evidence and the prequential Hyvarinen score (optional).
#'It is a wrapper of the function \code{smc_} with a time budgeting and partial save feature.
#'@export
smc = function(observations, model, algorithmic_parameters){
  # Set the time budget if needed
  saveprompt  = "not saved (no savefilename provided or option save is off)"
  if (!is.null(algorithmic_parameters$time_budget)){
    if (!is.null(algorithmic_parameters$save)&&(algorithmic_parameters$save == TRUE)&&!is.null(algorithmic_parameters$savefilename)){
      saveprompt  = paste("saved in",algorithmic_parameters$savefilename)
    }
    cat(strftime(Sys.time()),", time budget =",algorithmic_parameters$time_budget,"sec\n")
    setTimeLimit(elapsed = algorithmic_parameters$time_budget)
  }
  # Run the SMC with possible interruption
  results = tryCatch(smc_(observations, model, algorithmic_parameters),
                     error = function(e) {
                       if (regexpr("time limit",e$message) == -1) {print(e); return (NULL)}
                       else {cat("Time limit reached: partial results",saveprompt,"\n"); return (NULL)}
                     })
  # Resets time budget to infinity
  setTimeLimit()
  return (results)
}
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
smc_ = function(observations, model, algorithmic_parameters){
  # Parse algorithmic parameters and set flags accordingly
  nobservations = ncol(observations)
  Ntheta = algorithmic_parameters$Ntheta
  nmoves = algorithmic_parameters$nmoves
  resampling = algorithmic_parameters$resampling
  ess_objective = algorithmic_parameters$ess_threshold*algorithmic_parameters$Ntheta
  if (algorithmic_parameters$hscore) {observation_type = tolower(model$observation_type)}
  #-------------------------------------------------------------------------------------------------------
  # Monitor progress if needed
  if (algorithmic_parameters$progress) {
    print(paste("Started at:",Sys.time()))
    progbar = txtProgressBar(min = 0,max = nobservations,style=3)
    time_start = proc.time()
  }
  #-------------------------------------------------------------------------------------------------------
  # Initialize empty arrays and lists to store the results
  ESS = array(NA,dim = c(nobservations)) #ESS at successive times t
  incr_logevidence = array(NA,dim = c(nobservations)) #incremental log-evidence at successive times t
  incr_hscore = array(NA,dim = c(nobservations)) # OPTIONAL: incremental Hyvarinen score at successive times t
  incr_hscore_kde = array(NA,dim = c(nobservations)) # OPTIONAL: incremental Hyvarinen score at successive times t using kernel density estimators
  rejuvenation_times = c() #successive times where resampling is triggered
  rejuvenation_rate = c() #successive acceptance rates of resampling
  thetas_history = list() #successive sets of particles theta
  normw_history = list() #successive sets of normalized weights for theta
  byproducts_history = list() #successive byproducts (one for each theta at each time step)
  # # if we start from a proposal instead of the prior (e.g. improper prior)
  # # then the weights should be initialized differently
  # if (is.null(algorithmic_parameters$rinitial_theta)){
  #   thetas = model$rprior(Ntheta)
  #   logtargetdensities = model$dprior(thetas) # log target density at current particles
  #   normw = rep(1/Ntheta, Ntheta) # normalized weights
  #   logw = rep(0, Ntheta) # log normalized weights
  # } else {
  #   # this assumes that the posterior is proper after 1 observation.
  #   # For the general case where the posterior only becomes proper after k observations,
  #   # we should target directly the posterior at time k and start assimilating observations
  #   # from time (k+1) to T
  #   thetas = algorithmic_parameters$rinitial_theta(Ntheta)
  #   logtargetdensities = model$dprior(thetas) # log target density at current particles
  #   logw = logtargetdensities - algorithmic_parameters$dinitial_theta(thetas)
  #   w = exp(logw - max(logw))
  #   normw = w / sum(w)
  # }
  # Assuming the prior distribution is proper
  thetas = model$rprior(Ntheta)
  logtargetdensities = apply(thetas,2,model$dprior) # log target density at current particles
  normw = rep(1/Ntheta, Ntheta) # normalized weights
  logw = rep(0, Ntheta) # log normalized weights
  # save thetas if needed
  if (algorithmic_parameters$store_thetas_history){
    thetas_history[[1]] = thetas
    normw_history[[1]] = normw
  }
  #-------------------------------------------------------------------------------------------------------
  # initialize possible byproducts (e.g. Kalman filters, etc ...)
  byproducts = list()
  if (!is.null(model$initialize_byproducts)) {
    for (itheta in 1:Ntheta){
      byproducts[[itheta]] = model$initialize_byproducts(thetas[,itheta], observations)
    }
  } else {
    byproducts = NULL
  }
  #-------------------------------------------------------------------------------------------------------
  # Assimilate observations one by one
  for (t in 1:nobservations){
    #-------------------------------------------------------------------------------------------------------
    # OPTIONAL: compute the incremental hscore for discrete observations
    if (algorithmic_parameters$hscore && (observation_type=="discrete")) {
      incr_hscore[t] = hincrement_discrete_smc(thetas, normw, byproducts, t, observations, model,
                                               logtargetdensities, algorithmic_parameters)
    }
    #-------------------------------------------------------------------------------------------------------
    # OPTIONAL: compute the incremental hscore for continuous observations using kernel density estimators
    if (algorithmic_parameters$hscore && (observation_type=="continuous") && algorithmic_parameters$use_kde) {
      # compute incremental H score (with theta from time t-1)
      incr_hscore_kde[t] = hincrementContinuous_smc_kde(t, model, observations,thetas,normw,
                                                        byproducts, logtargetdensities,
                                                        algorithmic_parameters)
    }
    #-------------------------------------------------------------------------------------------------------
    # Assimilate the next observation
    results = assimilate_one_smc(thetas,byproducts,t,observations,model,logtargetdensities,logw,normw,algorithmic_parameters)
    #-------------------------------------------------------------------------------------------------------
    # Update the particles theta and compute the log-evidence
    thetas = results$thetas
    normw = results$normw
    logw = results$logw
    logtargetdensities = results$logtargetdensities
    incr_logevidence[t] = results$logcst
    if (!is.null(byproducts)) {byproducts = results$byproducts}
    #-------------------------------------------------------------------------------------------------------
    # OPTIONAL: compute the incremental hscore for continuous observations
    if (algorithmic_parameters$hscore && (observation_type=="continuous")) {
      incr_hscore[t] = hincrement_continuous_smc(thetas, normw, byproducts, t, observations, model,
                                                 logtargetdensities, algorithmic_parameters)
    }
    #-------------------------------------------------------------------------------------------------------
    # do some book-keeping
    ESS[t] = results$ESS
    if (!is.na(results$rejuvenation_time)) {
      rejuvenation_times = c(rejuvenation_times, results$rejuvenation_time)
    }
    if (!is.na(results$rejuvenation_rate)) {
      rejuvenation_rate = c(rejuvenation_rate, results$rejuvenation_rate)
    }
    if (algorithmic_parameters$store_thetas_history){
      thetas_history[[t+1]] = thetas
      normw_history[[t+1]] = normw
    }
    if (algorithmic_parameters$store_byproducts_history){
      byproducts_history[[t+1]] = byproducts
    }
    #-------------------------------------------------------------------------------------------------------
    # Update progress bar if needed
    if (algorithmic_parameters$progress) {
      setTxtProgressBar(progbar, t)
    }
    #-------------------------------------------------------------------------------------------------------
    # save partial results if needed
    if (algorithmic_parameters$save) {
      # save the variables required to resume and proceed further, in case of interrupted run
      # NOTE: thetas, normw, and PFs, are retrievable from the history stored in results_so_far
      required_to_resume = list(t = t, logw = logw, logtargetdensities = logtargetdensities,
                                observations = observations, model = model,
                                algorithmic_parameters = algorithmic_parameters)
      # save the results obtained up to this time
      results_so_far = list(thetas_history = thetas_history, normw_history = normw_history,
                            incr_logevidence = incr_logevidence[1:t], incr_hscore = incr_hscore[1:t],  ESS = ESS[1:t],
                            rejuvenation_times = rejuvenation_times, rejuvenation_rate = rejuvenation_rate,
                            method = 'SMC', incr_hscore_kde = incr_hscore_kde)
      # if the history of theta-particles is not saved, just keep the most recent ones
      if (!algorithmic_parameters$store_thetas_history){
        required_to_resume$thetas = thetas; required_to_resume$normw = normw
      }
      # if the history of byproducts is not saved, just keep the most recent ones
      if (algorithmic_parameters$store_byproducts_history){
        results_so_far$byproducts_history = byproducts_history
      }
      else {
        required_to_resume$byproducts = byproducts
      }
      # save into RDS file
      saveRDS(c(required_to_resume,results_so_far),file = algorithmic_parameters$savefilename)
    }
    #-------------------------------------------------------------------------------------------------------
  }
  # Update progress bar if needed
  if (algorithmic_parameters$progress) {
    close(progbar)
    time_end = proc.time()-time_start
    cat(paste("SMC: T = ",toString(nobservations),", Ntheta = ",toString(Ntheta),"\n",sep=""))
    print(time_end)
  }
  # If no need to store the latest particles or byproducts, set them to NULL before returning the results
  if (!algorithmic_parameters$store_last_thetas) {thetas = NULL; normw = NULL}
  if (!algorithmic_parameters$store_last_byproducts) {byproducts = NULL}
  # Return the results as a list
  return (list(thetas = thetas, normw = normw, byproducts = byproducts, logtargetdensities = logtargetdensities,
               thetas_history = thetas_history, normw_history = normw_history, byproducts_history = byproducts_history,
               logevidence = cumsum(incr_logevidence), hscore = cumsum(incr_hscore), hscoreKDE = cumsum(incr_hscore_kde),
               ESS = ESS, rejuvenation_times = rejuvenation_times, rejuvenation_rate = rejuvenation_rate,
               method = 'SMC', algorithmic_parameters = algorithmic_parameters))
}

