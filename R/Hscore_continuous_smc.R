# This function computes successive prequential Hyvarinen score for continuous observations by running SMC.
# It also computes the successive log-evidence as a by-product.
hscore_continuous_smc <- function(observations, model, algorithmic_parameters){
  # Set default values for the missing fields
  algorithmic_parameters = set_default_algorithmic_parameters(algorithmic_parameters)
  model = set_default_model(model)
  # Parse algorithmic parameters and set flags accordingly
  nobservations = ncol(observations)
  Ntheta = algorithmic_parameters$Ntheta
  nmoves = algorithmic_parameters$nmoves
  resampling = algorithmic_parameters$resampling
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
  # assimilate observations one by one
  for (t in 1:nobservations){
    results <- assimilate_one_smc(thetas, t, observations, model, Ntheta, ess_objective,
                                  nmoves, resampling, logtargetdensities, logw, normw,
                                  byproducts, algorithmic_parameters$verbose)
    thetas = results$thetas
    normw = results$normw
    logw = results$logw
    logtargetdensities = results$logtargetdensities
    logevidence[t] = results$logcst
    if (!is.null(byproducts)) {byproducts = results$byproducts}
    # compute prequential H score
    Hscore[t] = hincrementContinuous_smc(t, model, observations,thetas,normw,byproducts,Ntheta)
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
    cat(paste("Hscore: T = ",toString(nobservations),", Ntheta = ",toString(Ntheta),", Nx = ",toString(Nx),"\n",sep = ""))
    print(time_end)
  }
  return(list(thetas_history = thetas_history, normw_history = normw_history, logevidence = cumsum(logevidence),
              logtargetdensities = logtargetdensities, hscore = cumsum(Hscore),
              rejuvenation_times = rejuvenation_times[!is.na(rejuvenation_times)],
              rejuvenation_accept_rate = rejuvenation_accept_rate[!is.na(rejuvenation_accept_rate)]))
}
