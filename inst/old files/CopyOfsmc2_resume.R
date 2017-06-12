#'@rdname smc2_resume
#'@title smc2_resume
#'@description This function resumes a previously interrupted SMC2 run.
#'It is a wrapper of the function \code{smc2_resume_} with a time budgeting and partial save feature.
#'@export
smc2_resume = function(RDSsave=NULL, savefilename=NULL, next_observations=NULL, new_algorithmic_parameters=NULL){
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
  results = tryCatch(smc2_resume_(RDSsave, algorithmic_parameters, next_observations),
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
smc2_resume_ = function(RDSsave, algorithmic_parameters, next_observations=NULL){
  # Update the observations to assimilate
  observations = cbind(RDSsave$observations, next_observations)
  nobservations = ncol(observations)
  n_assimilated = RDSsave$t
  # Get the number of particles
  Ntheta = RDSsave$algorithmic_parameters$Ntheta
  # Retrieve the theta-particles history (possibly empty lists)
  thetas_history = RDSsave$thetas_history
  normw_history = RDSsave$normw_history
  # Retrieve the X-particles history by reconstructings trees
  if (algorithmic_parameters$store_X_history) {
    PF_history = lapply(1:n_assimilated,function(t)lapply(1:Ntheta,function(i)c(RDSsave$PF_history_no_tree[[t+1]][[i]], tree = tree_reconstruct(RDSsave$trees_attributes_history[[t+1]][[i]]))))
    PF_history = c(list(NULL),PF_history)
  } else {
    PF_history = list()
  }
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
  increase_Nx_times = RDSsave$increase_Nx_times
  increase_Nx_values = RDSsave$increase_Nx_values
  if (n_assimilated == nobservations){
    # all observations have been assimilated, we can return the result
    return(list(thetas_history = thetas_history, normw_history = normw_history, logtargetdensities = logtargetdensities,
                PF_history = PF_history, logevidence = cumsum(incr_logevidence), hscore = cumsum(incr_hscore),
                ESS = ESS, rejuvenation_times = rejuvenation_times, rejuvenation_rate = rejuvenation_rate,
                increase_Nx_times = increase_Nx_times, increase_Nx_values = increase_Nx_values,method = 'SMC2'))
  } else {
    # Get the latest particles and their normalized weights
    if (algorithmic_parameters$store_thetas_history) {
      thetas = thetas_history[[n_assimilated+1]]
      normw = normw_history[[n_assimilated+1]]
    } else {
      thetas = RDSsave$thetas
      normw = RDSsave$normw
    }
    # Get the latest particle filters
    if (algorithmic_parameters$store_X_history) {
      PFs = PF_history[[n_assimilated+1]]
    } else {
      PFs = lapply(1:Ntheta,function(i)c(RDSsave$PFs_no_trees[[i]],tree=tree_reconstruct(RDSsave$trees_attributes[[i]])))
    }
    # update remaining parameters and # Get all the variables required to resume SMC2 run
    Nx = PFs[[1]]$Nx
    logw = RDSsave$logw
    model = RDSsave$model
    if (algorithmic_parameters$hscore) {observation_type = tolower(model$observation_type)}
    #  OPTIONAL: For discrete hscore, retrieve array of most recent particles X targeting predictive distributions
    if (algorithmic_parameters$hscore && (observation_type=="discrete")) {
      Xprevious = RDSsave$Xprevious
      Xpred = RDSsave$Xpred
      XnormW_previous = RDSsave$XnormW_previous
    }
    # Monitor progress if needed
    if (algorithmic_parameters$progress) {
      print(paste("Started at:",Sys.time()))
      progbar = txtProgressBar(min = 0,max = nobservations,initial = n_assimilated,style=3)
      time_start = proc.time()
    }
    ########################################################################
    ################ Assimilate the remaining observations #################
    ########################################################################
    for (t in (n_assimilated+1):nobservations){
      # OPTIONAL: compute the incremental hscore for discrete observations
      if (algorithmic_parameters$hscore && (observation_type=="discrete")) {
        # Construct particles targeting the one-step-ahead predictive (need to reconstruct since size Nx might change)
        if (t > 1){
          Nx = PFs[[1]]$Nx
          Xpred = array(NA,dim = c(model$dimX,Nx,Ntheta)) # (need to reconstruct since size Nx might change)
          for (itheta in 1:Ntheta){
            X = Xprevious[,,itheta]
            if (is.null(dim(X))){
              Xpred[,,itheta] = model$rtransition(matrix(X,ncol = Nx), t, thetas[,itheta])
            }
            else{
              Xpred[,,itheta] = model$rtransition(X, t, thetas[,itheta])
            }
          }
        }
        # compute incremental H score (with theta from time t-1, see formula in the paper)
        incr_hscore[t] = Hd_smc2(t,model,observations[,t],thetas,normw,Xpred,XnormW_previous)
      }
      # Assimilate the next observation
      results = assimilate_one_smc2(thetas,PFs,t,observations,model,logtargetdensities,logw,normw,algorithmic_parameters)
      # Update the particles theta and compute the log-evidence
      thetas = results$thetas
      normw = results$normw
      logw = results$logw
      PFs = results$PFs
      logtargetdensities = results$logtargetdensities
      incr_logevidence[t] = results$logcst
      # OPTIONAL: compute incremental hscore here for continuous observations and update particles for discrete case
      if (algorithmic_parameters$hscore) {
        if (observation_type=="continuous") {
          incr_hscore[t] = hincrementContinuous_smc2(t,model,observations[,t,drop=FALSE],thetas,normw,PFs,Ntheta)
        } else if (observation_type=="discrete") {
          #matrix of normalized weights for X at previous step (need to reconstruct since size Nx might have changed)
          Nx = PFs[[1]]$Nx
          Xprevious = array(NA,dim = c(Nx, model$dimX, Ntheta))
          XnormW_previous = matrix(NA, nrow = PFs[[1]]$Nx, ncol = Ntheta)
          for (itheta in 1:Ntheta){
            Xprevious[,,itheta] = PFs[[itheta]]$X
            XnormW_previous[,itheta] = PFs[[itheta]]$xnormW
          }
        }
      }
      # do some book-keeping
      ESS[t] = results$ESS
      if (!is.na(results$rejuvenation_time)) {rejuvenation_times = c(rejuvenation_times, results$rejuvenation_time)}
      if (!is.na(results$rejuvenation_rate)) {rejuvenation_rate = c(rejuvenation_rate, results$rejuvenation_rate)}
      if (!is.na(results$increase_Nx_times)) {increase_Nx_times = c(increase_Nx_times, results$increase_Nx_times)}
      if (!is.na(results$increase_Nx_values)) {increase_Nx_values = c(increase_Nx_values, results$increase_Nx_values)}
      if (algorithmic_parameters$store_thetas_history){
        thetas_history[[t+1]] = thetas
        normw_history[[t+1]] = normw
      }
      if (algorithmic_parameters$store_X_history){
        PF_history[[t+1]] = PFs
      }
      # Update progress bar if needed
      if (algorithmic_parameters$progress) {
        setTxtProgressBar(progbar, t)
      }
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
        if (algorithmic_parameters$hscore && (observation_type=="discrete")) {
          # additional variables required for the discrete case
          required_for_discrete = list(Xpred = Xpred, XnormW_previous = XnormW_previous)
          saveRDS(c(required_to_resume,required_for_discrete,results_so_far),file = algorithmic_parameters$savefilename)
        } else {
          saveRDS(c(required_to_resume,results_so_far),file = algorithmic_parameters$savefilename)
        }
      }
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
                thetas_history = thetas_history, normw_history = normw_history, PF_history = PF_history,
                logevidence = cumsum(incr_logevidence), hscore = cumsum(incr_hscore),
                ESS = ESS, rejuvenation_times = rejuvenation_times, rejuvenation_rate = rejuvenation_rate,
                increase_Nx_times = increase_Nx_times, increase_Nx_values = increase_Nx_values,
                method = 'SMC2', algorithmic_parameters = algorithmic_parameters))
  }
}
