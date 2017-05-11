#'@rdname smc_sampler
#'@title smc_sampler
#'@description This function runs the SMC algorithm, using tempering.
#'@export
smc_sampler <- function(observations, model, algorithmic_parameters){
  Ntheta <- algorithmic_parameters$Ntheta
  nmoves = algorithmic_parameters$nmoves
  resampling = algorithmic_parameters$resampling
  ess_objective <- algorithmic_parameters$ess_threshold*algorithmic_parameters$Ntheta
  nobservations <- ncol(observations)
  # Set flags for progress tracking
  if (is.null(algorithmic_parameters$progress)) algorithmic_parameters$progress = FALSE
  if (is.null(algorithmic_parameters$verbose)) algorithmic_parameters$verbose = FALSE
  if (algorithmic_parameters$progress) {
    print(paste("Started at:",Sys.time()))
    progbar = txtProgressBar(min = 0,max = nobservations,style=3)
    count = 0
    time_start = proc.time()
  }
  ESS = array(NA,dim = c(nobservations)) #ESS at successive times t
  logevidence = array(NA,dim = c(nobservations)) #log-evidence at successive times t
  rejuvenation_times <- c() #successive times where resampling is triggered
  rejuvenation_accept_rate <- c() #successive acceptance rates of resampling
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
  #
  thetas_history <- list()
  thetas_history[[1]] <- thetas
  normw_history <- list()
  normw_history[[1]] <- normw
  #
  logevidence <- rep(0, nobservations)
  #
  for (t in 1:nobservations){
    results <- assimilate_one_smc(thetas, t, observations, model, Ntheta, ess_objective,
                                  nmoves, resampling, logtargetdensities, logw, normw,
                                  algorithmic_parameters$verbose)
    thetas <- results$thetas
    normw <- results$normw
    logw <- results$logw
    logtargetdensities <- results$logtargetdensities
    logevidence[t] <- results$logcst
    thetas_history[[t+1]] <- thetas
    normw_history[[t+1]] <- normw
    # Update progress bar if needed
    if (algorithmic_parameters$progress) {
      count = count + 1
      setTxtProgressBar(progbar, count)
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
              logevidence = cumsum(logevidence), logtargetdensities = logtargetdensities))
}

