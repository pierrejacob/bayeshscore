#'@rdname smc2
#'@title smc2
#'@description This function runs the SMC2 algorithm, using adapive Nx and tempering.
#'It also computes the log-evidence and the prequential Hyvarinen score (optional).
#'It is a wrapper of the function \code{smc2_} with a time budgeting and partial save feature.
#'@export
smc2 = function(observations, model, algorithmic_parameters){
  # load TreeClass if needed
  tryCatch(TreeClass, error = function(e) {
    if (regexpr("TreeClass",e$message) > -1) {module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM");
    TreeClass <<- module_tree$Tree}})
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
  results = tryCatch(smc2_(observations, model, algorithmic_parameters),
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
smc2_ = function(observations, model, algorithmic_parameters){
  # Parse algorithmic parameters and set flags accordingly
  nobservations = ncol(observations)
  Ntheta = algorithmic_parameters$Ntheta
  Nx = algorithmic_parameters$Nx
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
  rejuvenation_times = c() #successive times where resampling is triggered
  rejuvenation_rate = c() #successive acceptance rates of resampling
  increase_Nx_times = c() #successive times where increasing Nx is triggered
  increase_Nx_values = c() #successive values of Nx
  thetas_history = list() #successive sets of particles theta
  normw_history = list() #successive sets of normalized weights for theta
  logtargetdensities_history = list() #successive target log-densities for each particle theta
  PF_history = list() #successive particle filters (one for each theta at each time step)
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
  if (algorithmic_parameters$store_thetas_history){
    thetas_history[[1]] = thetas
    normw_history[[1]] = normw
    logtargetdensities_history[[1]] = logtargetdensities
  }
  #-------------------------------------------------------------------------------------------------------
  # Initialize filters (first observation passed as argument just to initialize the fields of PF)
  PFs = list() # list of particle filters (one for each theta)
  for (itheta in 1:Ntheta){
    theta = thetas[,itheta]
    PFs[[itheta]] = conditional_particle_filter(observations[,1,drop = FALSE], model, theta, Nx)
    #Note: the CPF performs a regular PF when no conditioning path is provided
  }
  #-------------------------------------------------------------------------------------------------------
  # Assimilate observations one by one
  for (t in 1:nobservations){
    #-------------------------------------------------------------------------------------------------------
    # OPTIONAL: compute the incremental hscore for discrete observations
    if (algorithmic_parameters$hscore && (observation_type=="discrete")) {
      incr_hscore[t] = hincrement_discrete_smc2(thetas, normw, PFs, t, observations, model, logtargetdensities, algorithmic_parameters)
    }
    #-------------------------------------------------------------------------------------------------------
    # Assimilate the next observation
    results = assimilate_one_smc2(thetas,PFs,t,observations,model,logtargetdensities,logw,normw,algorithmic_parameters)
    #-------------------------------------------------------------------------------------------------------
    # Update the particles theta and compute the log-evidence
    thetas = results$thetas
    normw = results$normw
    logw = results$logw
    PFs = results$PFs
    logtargetdensities = results$logtargetdensities
    incr_logevidence[t] = results$logcst
    #-------------------------------------------------------------------------------------------------------
    # OPTIONAL: compute incremental hscore here for continuous observations and update particles for discrete case
    if (algorithmic_parameters$hscore && observation_type=="continuous") {
      incr_hscore[t] = hincrement_continuous_smc2(thetas, normw, PFs, t, observations, model, logtargetdensities, algorithmic_parameters)
    }
    #-------------------------------------------------------------------------------------------------------
    # do some book-keeping
    ESS[t] = results$ESS
    if (!is.na(results$rejuvenation_time)) {rejuvenation_times = c(rejuvenation_times, results$rejuvenation_time)}
    if (!is.na(results$rejuvenation_rate)) {rejuvenation_rate = c(rejuvenation_rate, results$rejuvenation_rate)}
    if (!is.na(results$increase_Nx_times)) {increase_Nx_times = c(increase_Nx_times, results$increase_Nx_times)}
    if (!is.na(results$increase_Nx_values)) {increase_Nx_values = c(increase_Nx_values, results$increase_Nx_values)}
    if (algorithmic_parameters$store_thetas_history){
      thetas_history[[t+1]] = thetas
      normw_history[[t+1]] = normw
      logtargetdensities_history[[t+1]] = logtargetdensities
    }
    if (algorithmic_parameters$store_X_history){
      PF_history[[t+1]] = PFs
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
                            logtargetdensities_history = logtargetdensities_history,
                            incr_logevidence = incr_logevidence[1:t], incr_hscore = incr_hscore[1:t],  ESS = ESS[1:t],
                            rejuvenation_times = rejuvenation_times, rejuvenation_rate = rejuvenation_rate,
                            increase_Nx_times = increase_Nx_times, increase_Nx_values = increase_Nx_values,
                            method = 'SMC2')
      # if the history of x-particles is not saved, just keep the most recent ones
      if (algorithmic_parameters$store_X_history){
        results_so_far$PF_history_no_tree = lapply(1:(t+1),function(j)lapply(1:Ntheta,function(i)PF_history[[j]][[i]][names(PF_history[[j]][[i]])!="tree"]))
        results_so_far$trees_attributes_history = lapply(1:(t+1),function(j)trees_getattributes(lapply(1:Ntheta,function(i)PF_history[[j]][[i]]$tree)))
      } else {
        required_to_resume$PFs_no_trees = lapply(1:Ntheta,function(i)PFs[[i]][names(PFs[[i]])!="tree"])
        required_to_resume$trees_attributes = trees_getattributes(lapply(1:Ntheta,function(i)PFs[[i]]$tree))
      }
      # if the history of theta-particles is not saved, just keep the most recent ones
      if (!algorithmic_parameters$store_thetas_history){
        required_to_resume$thetas = thetas
        required_to_resume$normw = normw
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
    cat(paste("SMC2: T = ",toString(nobservations),", Ntheta = ",toString(Ntheta),", Nx (last) = ",toString(PFs[[1]]$Nx),"\n",sep=""))
    print(time_end)
  }
  # If no need to store the latest particles or byproducts, set them to NULL before returning the results
  if (!algorithmic_parameters$store_last_thetas) {thetas = NULL; normw = NULL}
  if (!algorithmic_parameters$store_last_X) {PFs = NULL}
  # Return the results as a list
  return(list(thetas = thetas, normw = normw, PFs = PFs, logtargetdensities = logtargetdensities,
              thetas_history = thetas_history, normw_history = normw_history,
              logtargetdensities_history = logtargetdensities_history, PF_history = PF_history,
              logevidence = cumsum(incr_logevidence), hscore = cumsum(incr_hscore),
              ESS = ESS, rejuvenation_times = rejuvenation_times, rejuvenation_rate = rejuvenation_rate,
              increase_Nx_times = increase_Nx_times, increase_Nx_values = increase_Nx_values,
              method = 'SMC2', algorithmic_parameters = algorithmic_parameters))
}
