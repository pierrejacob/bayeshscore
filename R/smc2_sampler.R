#'@rdname smc2_sampler
#'@title smc2_sampler
#'@description This function runs the SMC^2 algorithm, using adapive Nx and tempering.
#'@export
smc2_sampler <- function(observations, model, algorithmic_parameters){
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
    progbar = txtProgressBar(min = 0,max = nobservations,style=3)
    time_start = proc.time()
  }
  # Initialize empty arrays and lists to store the results
  ESS = array(NA,dim = c(nobservations)) #ESS at successive times t
  Hscore = array(NA,dim = c(nobservations)) #prequential Hyvarinen score at successive times t
  logevidence = array(NA,dim = c(nobservations)) #log-evidence at successive times t
  thetas_history = list() #successive sets of particles theta
  normw_history = list() #successive sets of normalized weights for theta
  PF_history = list() #successive particle filters (one for each theta at each time step)
  # # Initialize SMC2 by sampling from prior (default), or from specified proposal
  # # (e.g. relevant in case the prior is vague)
  # if (is.null(algorithmic_parameters$rinitial_theta)){
  thetas = model$rprior(Ntheta)
  # } else {
  # thetas = algorithmic_parameters$rinitial_theta(Ntheta)
  # }
  logtargetdensities = apply(thetas, 2, model$dprior) # log target density evaluations at current particles
  normw = rep(1/Ntheta, Ntheta) # normalized weights
  logw = rep(0, Ntheta) # log normalized weights
  PFs = list() # list of particle filters (one for each theta)
  thetas_history[[1]] = thetas
  normw_history[[1]] = normw
  # Initialize filters (first observation passed as argument just to initialize the fields of PF)
  for (itheta in 1:Ntheta){
    theta = thetas[,itheta]
    PFs[[itheta]] = conditional_particle_filter(observations[,1,drop = FALSE], model, theta, Nx)
    #the CPF performs a regular PF when no conditioning path is provided
  }
  # Assimilate observations one by one
  for (t in 1:nobservations){
    results = assimilate_one_smc2(thetas, PFs, t, observations, model, Ntheta, ess_objective,
                              nmoves, resampling, logtargetdensities, logw, normw,
                              algorithmic_parameters$verbose)
    thetas = results$thetas
    normw = results$normw
    logw = results$logw
    PFs = results$PFs
    logtargetdensities = results$logtargetdensities
    logevidence[t] = results$logcst
    thetas_history[[t+1]] = thetas
    normw_history[[t+1]] = normw
    if (algorithmic_parameters$store_X){
      PF_history[[t+1]] = PFs
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
    cat(paste("T = ",toString(nobservations),", Ntheta = ",toString(Ntheta),", Nx = ",toString(Nx),"\n",sep=""))
    print(time_end)
  }
  return(list(thetas_history = thetas_history, normw_history = normw_history, PF_history = PF_history,
              logevidence = cumsum(logevidence), logtargetdensities = logtargetdensities))
}
