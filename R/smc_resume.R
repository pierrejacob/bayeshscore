#'@rdname smc_resume
#'@title smc_resume
#'@description This function resumes a previously interrupted SMC run.
#'It is a wrapper of the function \code{smc_resume_} with a time budgeting and partial save feature.
#'@export
smc_resume = function(RDSsave=NULL, savefilename=NULL, next_observations=NULL, new_algorithmic_parameters=NULL){
  ########################################################################
  ################ Load partial results from save file ###################
  ########################################################################
  # load RDS file
  if (is.null(RDSsave)) {
    if (is.null(savefilename)) {
      cat("partial results (loaded from RDS file) or path to RDS file (savefilename) must be provided\n")
      return (NULL)
    } else {
      RDSsave = readRDS(savefilename)
    }
  }
  # update new algorithmic parameters and flags. NOTE: some parameters CANNOT be modified (e.g. Ntheta, model, ...)
  algorithmic_parameters = RDSsave$algorithmic_parameters
  mutable = c("progress","verbose","save","savefilename","time_budget","ess_threshold","nmoves","resampling","proposalmove")
  for (i in 1:length(mutable)) {
    if (!is.null(new_algorithmic_parameters[[mutable[i]]])) {
      algorithmic_parameters[[mutable[i]]] = new_algorithmic_parameters[[mutable[i]]]
    }
  }
  # Set the time budget if needed
  if (!is.null(algorithmic_parameters$time_budget)){
    if (is.null(algorithmic_parameters$save)||(algorithmic_parameters$save == FALSE)||is.null(algorithmic_parameters$savefilename)){
      saveprompt  = "not saved (no savefilename provided or option save is off)"
    } else {
      saveprompt  = paste("saved in",algorithmic_parameters$savefilename)
    }
    cat(strftime(Sys.time()),", time budget =",algorithmic_parameters$time_budget,"sec\n")
    setTimeLimit(elapsed = algorithmic_parameters$time_budget)
  }
  # Run the SMC with possible interruption
  results = tryCatch(smc_resume_(RDSsave, algorithmic_parameters, next_observations),
                     error = function(e) {cat("Time limit reached: partial results",saveprompt,"\n"); NULL})
  # Resets time budget to infinity
  setTimeLimit()
  return (results)
}
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
smc_resume_ = function(RDSsave, algorithmic_parameters, next_observations=NULL){
  # Update the observations to assimilate
  observations = cbind(RDSsave$observations, next_observations)
  nobservations = ncol(observations)
  n_assimilated = RDSsave$t
  # Get the number of particles
  Ntheta = RDSsave$algorithmic_parameters$Ntheta
  # Retrieve the theta-particles history and byproducts history(possibly empty lists)
  thetas_history = RDSsave$thetas_history
  normw_history = RDSsave$normw_history
  byproducts_history = RDSsave$byproducts_history
  # Retrieve the results up to now
  logtargetdensities = RDSsave$logtargetdensities # log target density evaluations at latest particles
  ESS = array(NA,dim = c(nobservations)) #ESS at successive times
  ESS[1:n_assimilated] = RDSsave$ESS #reload previous ESS
  incr_logevidence = array(NA,dim = c(nobservations)) #incremental log-evidence at successive times
  incr_logevidence[1:n_assimilated] = RDSsave$incr_logevidence #reload previous logevidence
  incr_hscore = array(NA,dim = c(nobservations))
  if (algorithmic_parameters$hscore) {
    # OPTIONAL: incremental Hyvarinen score
    incr_hscore[1:n_assimilated] = RDSsave$incr_hscore
  }
  # Retrieve the diagnostics up to now
  rejuvenation_times = RDSsave$rejuvenation_times
  rejuvenation_rate = RDSsave$rejuvenation_rate
  if (n_assimilated == nobservations){
    # all observations have been assimilated, we can return the result
    return (list(thetas_history = thetas_history, normw_history = normw_history, logtargetdensities = logtargetdensities,
                 byproducts_history = byproducts_history, logevidence = cumsum(incr_logevidence), hscore = cumsum(incr_hscore),
                 ESS = ESS, rejuvenation_times = rejuvenation_times, rejuvenation_rate = rejuvenation_rate))
  }
  # Get the latest particles and their normalized weights
  if (algorithmic_parameters$store_theta) {
    thetas = thetas_history[[n_assimilated+1]]
    normw = normw_history[[n_assimilated+1]]
  } else {
    thetas = RDSsave$thetas
    normw = RDSsave$normw
  }
  # Get the latest byproducts
  if (algorithmic_parameters$store_byproducts) {
    byproducts = RDSsave$byproducts_history[[n_assimilated+1]]
  } else {
    byproducts = RDSsave$byproducts
  }
  # update remaining parameters and # Get all the variables required to resume SMC run
  logw = RDSsave$logw
  model = RDSsave$model
  if (algorithmic_parameters$hscore) {observation_type = tolower(model$observation_type)}
  # Monitor progress if needed
  if (algorithmic_parameters$progress) {
    print(paste("Started at:",Sys.time()))
    progbar = txtProgressBar(min = 0,max = nobservations,initial = n_assimilated,style=3)
    time_start = proc.time()
  }
  #-------------------------------------------------------------------------------------------------------
  ########################################################################
  ################ Assimilate the remaining observations #################
  ########################################################################
  for (t in (n_assimilated+1):nobservations){
    #-------------------------------------------------------------------------------------------------------
    # OPTIONAL: compute the incremental hscore for discrete observations
    if (algorithmic_parameters$hscore && (observation_type=="discrete")) {
      # compute incremental H score (with theta from time t-1, see formula in the paper)
      incr_hscore[t] = Hd_smc(t,model,observations,thetas,normw,Ntheta,byproducts)
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
      incr_hscore[t] = hincrementContinuous_smc(t, model, observations,thetas,normw,byproducts,Ntheta)
    }
    #-------------------------------------------------------------------------------------------------------
    # do some book-keeping
    ESS[t] = results$ESS
    if (!is.na(results$rejuvenation_time)) {rejuvenation_times = c(rejuvenation_times, results$rejuvenation_time)}
    if (!is.na(results$rejuvenation_rate)) {rejuvenation_rate = c(rejuvenation_rate, results$rejuvenation_rate)}
    if (algorithmic_parameters$store_theta){thetas_history[[t+1]] = thetas; normw_history[[t+1]] = normw}
    if (algorithmic_parameters$store_byproducts){byproducts_history[[t+1]] = byproducts}
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
                            rejuvenation_times = rejuvenation_times, rejuvenation_rate = rejuvenation_rate)
      # if the history of theta-particles is not saved, just keep the most recent ones
      if (!algorithmic_parameters$store_theta){
        required_to_resume$thetas = thetas; required_to_resume$normw = normw
      }
      # if the history of byproducts is not saved, just keep the most recent ones
      if (algorithmic_parameters$store_byproducts){
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
  return (list(thetas_history = thetas_history, normw_history = normw_history, logtargetdensities = logtargetdensities,
               byproducts_history = byproducts_history, logevidence = cumsum(incr_logevidence), hscore = cumsum(incr_hscore),
               ESS = ESS, rejuvenation_times = rejuvenation_times, rejuvenation_rate = rejuvenation_rate))
}













