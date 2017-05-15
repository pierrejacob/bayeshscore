#'@rdname smc_sampler
#'@title smc_sampler
#'@description This function runs the SMC algorithm, using tempering.
#'@export
smc_sampler <- function(observations, model, algorithmic_parameters){
  ## Set default values for the missing fields
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
    progbar = txtProgressBar(min = 0,max = nobservations,style=3)
    time_start = proc.time()
  }
  # Initialize empty arrays and lists to store the results
  ESS = array(NA,dim = c(nobservations)) #ESS at successive times t
  logevidence = array(NA,dim = c(nobservations)) #log-evidence at successive times t
  thetas_history = list() #successive sets of particles theta
  normw_history = list() #successive sets of normalized weights for theta
  # if (is.null(algorithmic_parameters$rinitial_theta)){
  thetas <- model$rprior(Ntheta)
  # } else {
  # thetas <- algorithmic_parameters$rinitial_theta(Ntheta)
  # }
  # log target density evaluations at current particles
  logtargetdensities <- apply(thetas, 1, model$dprior)
  # normalized weights
  normw <- rep(1/Ntheta, Ntheta)
  logw <- rep(0, Ntheta)

  ########## if we start from a proposal instead of the prior (e.g. improper prior)
  ########## then the weights should be initialized differently:
  ########## log(prior_density) - log(proposal_density) ???

  thetas_history[[1]] <- thetas
  normw_history[[1]] <- normw
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
    thetas <- results$thetas
    normw <- results$normw
    logw <- results$logw
    logtargetdensities <- results$logtargetdensities
    logevidence[t] <- results$logcst
    if (!is.null(byproducts)) {byproducts = results$byproducts}
    thetas_history[[t+1]] <- thetas
    normw_history[[t+1]] <- normw
    # Update progress bar if needed
    if (algorithmic_parameters$progress) {
      setTxtProgressBar(progbar, t)
    }
  }
  # Update progress bar if needed
  if (algorithmic_parameters$progress) {
    close(progbar)
    time_end = proc.time()-time_start
    cat(paste("T = ",toString(nobservations),", Ntheta = ",toString(Ntheta),"\n",sep=""))
    print(time_end)
  }
  return(list(thetas_history = thetas_history, normw_history = normw_history,
              logevidence = cumsum(logevidence), logtargetdensities = logtargetdensities,
              byproducts = byproducts))
}

