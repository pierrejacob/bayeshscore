#'@rdname PMMH
#'@title PMMH
#'@description This function uses PMMH to sample from the posterior p(theta|Y1:t)
#'@export
PMMH = function(observations, model, algorithmic_parameters){
  M = algorithmic_parameters$M
  Nx = algorithmic_parameters$Nx
  burnin = algorithmic_parameters$burnin
  rMHproposal = algorithmic_parameters$rMHproposal
  dMHproposal = algorithmic_parameters$dMHproposal
  dimtheta = model$dimtheta
  if (algorithmic_parameters$progress) {
    print(paste("Started at:",Sys.time()))
    progbar = txtProgressBar(min = 0,max = M,style=3)
    count = 0
    time_start = proc.time()
  }
  #initialize
  thetas = matrix(NA,nrow = M,ncol = dimtheta)
  thetas[1,] = model$rprior(1)
  smc = bootstrap_particle_filter(observations, model, thetas[1,], algorithmic_parameters)
  log_p_hat_old = smc$log_p_y_hat
  if (algorithmic_parameters$progress) {
    count = count + 1
    setTxtProgressBar(progbar, count)
  }
  #iterate MH
  for (i in 2:M) {
    theta_old = matrix(thetas[i-1,], ncol = dimtheta)
    theta_new = rMHproposal(theta_old)
    if (model$dprior(theta_new,log = TRUE) == -Inf){
      thetas[i,] = theta_old
      next
    }
    else {
      smc = bootstrap_particle_filter(observations, model, theta_new, algorithmic_parameters)
      log_p_hat_new = smc$log_p_y_hat
      log_acceptance_num = model$dprior(theta_new,log = TRUE) + log_p_hat_new + dMHproposal(theta_old,theta_new)
      log_acceptance_den = model$dprior(theta_old,log = TRUE) + log_p_hat_old + dMHproposal(theta_new,theta_old)
      log_acceptance = log_acceptance_num - log_acceptance_den
      log_u = log(runif(1,min = 0,max = 1))
      if (log_u < log_acceptance){
        thetas[i,] = theta_new
        log_p_hat_old = log_p_hat_new
      }
      else {
        thetas[i,] = theta_old
      }
    }
    if (algorithmic_parameters$progress) {
      count = count + 1
      setTxtProgressBar(progbar, count)
    }
  }
  if (algorithmic_parameters$progress) {
    close(progbar)
    time_end = proc.time()-time_start
    cat(paste("PMMH: ",", M = ",toString(M), ", Nx = ",toString(Nx),"\n",sep = ""))
    print(time_end)
  }
  return (thetas[(burnin+1):M,])
}
