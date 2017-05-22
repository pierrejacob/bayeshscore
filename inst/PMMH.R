#'@rdname PMMH
#'@title PMMH
#'@description This function uses PMMH to sample from the posterior p(theta|Y_1:t)
#'@export
PMMH = function(observations, model, algorithmic_parameters){
  M = algorithmic_parameters$M # number of iterations
  Nx = algorithmic_parameters$Nx # number of particles
  burnin = algorithmic_parameters$burnin
  rMHproposal = algorithmic_parameters$rMHproposal # proposal sampler
  dMHproposal = algorithmic_parameters$dMHproposal # proposal density function
  dimtheta = model$dimtheta
  if (algorithmic_parameters$progress) {
    print(paste("Started at:",Sys.time()))
    progbar = txtProgressBar(min = 0,max = M,style=3)
    count = 0
    time_start = proc.time()
  }
  #initialize
  thetas = matrix(NA,nrow = dimtheta,ncol = M)
  thetas[,1] = model$rprior(1)
  log_p_hat_old = bootstrap_particle_filter(observations, model, thetas[,1], algorithmic_parameters)$log_p_y_hat
  log_prior_old = model$dprior(thetas[,1],log = TRUE)
  if (algorithmic_parameters$progress) {
    setTxtProgressBar(progbar, 1)
  }
  #iterate MH
  for (i in 2:M) {
    theta_old = thetas[,i-1,drop=FALSE]
    theta_new = rMHproposal(theta_old)
    log_prior_new = model$dprior(theta_new,log = TRUE)
    if (log_prior_new == -Inf){
      thetas[,i] = theta_old
      setTxtProgressBar(progbar, i)
      next
    }
    else {
      log_p_hat_new = bootstrap_particle_filter(observations, model, theta_new, algorithmic_parameters)$log_p_y_hat
      log_acceptance_num = log_prior_new + log_p_hat_new + dMHproposal(theta_old,theta_new)
      log_acceptance_den = log_prior_old + log_p_hat_old + dMHproposal(theta_new,theta_old)
      if (log(runif(1)) < (log_acceptance_num - log_acceptance_den)){
        thetas[,i] = theta_new
        log_p_hat_old = log_p_hat_new
        log_prior_old = log_prior_new
      }
      else {
        thetas[,i] = theta_old
      }
    }
    if (algorithmic_parameters$progress) {
      setTxtProgressBar(progbar, i)
    }
  }
  if (algorithmic_parameters$progress) {
    close(progbar)
    time_end = proc.time()-time_start
    cat(paste("PMMH: ",", M = ",toString(M), ", Nx = ",toString(Nx),"\n",sep = ""))
    print(time_end)
  }
  return (thetas[,(burnin+1):M,drop=FALSE])
}
