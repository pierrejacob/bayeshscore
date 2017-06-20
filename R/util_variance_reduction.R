#------------------------------------------------------------------------------#
#-----------  Some functions to improve particle-based estimation  ------------#
#------------------------------------------------------------------------------#
#'@rdname get_additional_particles_smc2
#'@title get_additional_particles_smc2
#'@description (SMC2 version) Generates additional particles given current ones, in order to reduce the variance
#'of the associated particle estimators. Nc is the desired total number of particles
#'to be used in the computation / estimation (so about Nc - Ntheta new particles are generated,
#'used in the estimation, and then discarded). The generated particles are equally weigthed.
#'Note that in the continuous observations case, the particles (theta,xt) target the posterior-filtering
#'distribution p(theta, xt | y_1, ..., y_t), whereas in the discrete case, they must target
#'the posterior-one-step-predictive distribution p(theta, xt | y_1, ..., y_(t-1)).
#'@export
get_additional_particles_smc2 = function(Nc, thetas, normw, PFs, t, observations, model,
                                         logtargetdensities, algorithmic_parameters, Ncx = NULL) {
  Ntheta = ncol(thetas)
  # each moves generates Ntheta new particles thetas.
  # By convention, the first move corresponds to the current particles thetas, so that nmoves = 1
  # even when no additional moves are needed
  nmoves = ceiling(Nc/Ntheta)
  # keep all the Ntheta*nmoves particles
  thetas_larger_sample = matrix(NA, nrow = nrow(thetas), ncol = nmoves*Ntheta)
  PFs_larger_sample = vector("list",nmoves*Ntheta)
  # Resample and move
  proposalmove = algorithmic_parameters$proposalmove(thetas,normw,model)
  # Resample particles theta according to normalized weights
  resampled_index = algorithmic_parameters$resampling(normw)
  thetas = thetas[,resampled_index,drop=FALSE]
  PFs = PFs[resampled_index]
  logtargetdensities = logtargetdensities[resampled_index]
  # By convention, the first move corresponds to the current particles thetas.
  thetas_larger_sample[,1:Ntheta] = thetas
  # Generate more x-particles for the current thetas if needed
  if (!is.null(Ncx)) {
    if (Ncx > PFs[[1]]$Nx) {
      for (i in 1:Ntheta) {
        PFs[[i]] = conditional_particle_filter(observations[,1:t,drop=FALSE], model, thetas[,i], Ncx)
      }
    }
  }
  Nx = PFs[[1]]$Nx
  PFs_larger_sample[1:Ntheta] = PFs
  if (algorithmic_parameters$verbose) {
    Ntheta_total = nmoves*Ntheta
    Ntheta_current = Ntheta
    cat("\r Additional particles generated:",Ntheta_current,"out of",Ntheta_total,"(",floor(Ntheta_current/Ntheta_total*100),"%)")
    flush.console()
  }
  # Generate more theta-particles (and associated x-particles) if needed
  if (nmoves >= 2){
    for (imove in 2:nmoves){
      theta_new_all = proposalmove$r(Ntheta)
      log_proposal_density_new_all = proposalmove$d(theta_new_all,log=TRUE)
      log_proposal_density_current = proposalmove$d(thetas,log=TRUE)
      accepts = 0
      for (i in 1:Ntheta) {
        theta_new = theta_new_all[,i]
        logprior_theta_new = model$dprior(theta_new, log = TRUE)
        if (is.infinite(logprior_theta_new)){
          thetas_larger_sample[,((imove-1)*Ntheta+i)] = thetas[,i]
          PFs_larger_sample[[((imove-1)*Ntheta+i)]] = PFs[[i]]
          next
        } else {
          #the CPF function performs a regular PF when no conditioning path is provided
          PF_new = conditional_particle_filter(observations[,1:t,drop=FALSE], model, theta_new, Nx)
          incremental_ll_new = PF_new$incremental_ll
          loglikelihood_new = incremental_ll_new[t]
          if (t > 1){
            loglikelihood_new = loglikelihood_new + sum(incremental_ll_new[1:(t-1)])
          }
          logtarget_new = logprior_theta_new + loglikelihood_new
          lognum = logtarget_new + log_proposal_density_current[i]
          logdenom = logtargetdensities[i] + log_proposal_density_new_all[i]
          logacceptance = lognum - logdenom
          logu = log(runif(1))
          if (logu <= logacceptance){
            thetas[,i] = theta_new
            logtargetdensities[i] = logtarget_new
            PFs[[i]] = PF_new
          }
          # otherwise do nothing (i.e. keep the current particles)
        }
        thetas_larger_sample[,((imove-1)*Ntheta+i)] = thetas[,i]
        PFs_larger_sample[[((imove-1)*Ntheta+i)]] = PFs[[i]]
      }
      if (algorithmic_parameters$verbose) {
        Ntheta_current = Ntheta_current + Ntheta
        cat("\r Additional particles generated:",Ntheta_current,"out of",Ntheta_total,"(",floor(Ntheta_current/Ntheta_total*100),"%)")
        flush.console()
      }
    }
  }
  if (algorithmic_parameters$verbose) {
    cat("\n Done. \n")
  }
  return (list(thetas = thetas_larger_sample, PFs = PFs_larger_sample))
}


#'@rdname get_additional_particles_smc
#'@title get_additional_particles_smc
#'@description (SMC version) Generates additional particles given current ones, in order to reduce the variance
#'of the associated particle estimators. Nc is the desired total number of particles
#'to be used in the computation / estimation (so about Nc - Ntheta new particles are generated,
#'used in the estimation, and then discarded). The generated particles are equally weigthed.
#'@export
get_additional_particles_smc = function(Nc, thetas, normw, byproducts, t, observations, model,
                                        logtargetdensities, algorithmic_parameters) {
  Ntheta = ncol(thetas)
  # each moves generates Ntheta new particles thetas.
  # By convention, the first move corresponds to the current particles thetas, so that nmoves = 1
  # even when no additional moves are needed
  nmoves = ceiling(Nc/Ntheta)
  # keep all the Ntheta*nmoves particles
  thetas_larger_sample = matrix(NA, nrow = nrow(thetas), ncol = nmoves*Ntheta)
  byproducts_larger_sample = vector("list",nmoves*Ntheta)
  # Resample and move
  proposalmove = algorithmic_parameters$proposalmove(thetas,normw,model)
  # Resample particles theta according to normalized weights
  resampled_index = algorithmic_parameters$resampling(normw)
  thetas = thetas[,resampled_index,drop=FALSE]
  if (!is.null(byproducts)) {byproducts = byproducts[resampled_index]}
  logtargetdensities = logtargetdensities[resampled_index]
  # By convention, the first move corresponds to the current particles thetas.
  thetas_larger_sample[,1:Ntheta] = thetas
  if (!is.null(byproducts)) {byproducts_larger_sample[1:Ntheta] = byproducts}
  if (algorithmic_parameters$verbose) {
    Ntheta_total = nmoves*Ntheta
    Ntheta_current = Ntheta
    cat("\r Additional particles generated:",Ntheta_current,"out of",Ntheta_total,"(",floor(Ntheta_current/Ntheta_total*100),"%)")
    flush.console()
  }
  if (nmoves >= 2){
    for (imove in 2:nmoves){
      thetas_new_all = proposalmove$r(Ntheta)
      log_proposal_density_new_all = proposalmove$d(thetas_new_all,log=TRUE)
      log_proposal_density_current = proposalmove$d(thetas,log=TRUE)
      accepts = 0
      for (i in 1:Ntheta) {
        theta_new = thetas_new_all[,i]
        logprior_theta_new = model$dprior(theta_new, log = TRUE)
        if (is.infinite(logprior_theta_new)){
          thetas_larger_sample[,((imove-1)*Ntheta+i)] = thetas[,i]
          if (!is.null(byproducts)) {byproducts_larger_sample[[((imove-1)*Ntheta+i)]] = byproducts[[i]]}
          next
        } else {
          # compute the log-likelihood of proposed theta using byproduct (e.g. Kalman filter) or analytically via the model
          if (!is.null(byproducts)){
            byproduct = model$initialize_byproducts(theta_new, observations)
            for (j in 1:t){
              byproduct = model$update_byproduct(byproduct, j, theta_new, observations)
            }
            incremental_ll_new = model$dpredictive(observations,t,theta_new,byproduct,log = TRUE)
          } else {
            incremental_ll_new = model$dpredictive(observations,t,theta_new,log = TRUE)
          }
          loglikelihood_new = incremental_ll_new
          if (t > 1){
            # compute the log-likelihood of proposed theta using byproduct (e.g. Kalman filter) or analytically via the model
            if (!is.null(byproducts)){
              loglikelihood_new = loglikelihood_new + model$likelihood(observations,t-1,theta_new,byproduct,log = TRUE)
            } else {
              loglikelihood_new = loglikelihood_new + model$likelihood(observations,t-1,theta_new,log = TRUE)
            }
          }
          logtarget_new = logprior_theta_new + loglikelihood_new
          lognum = logtarget_new + log_proposal_density_current[i]
          logdenom = logtargetdensities[i] + log_proposal_density_new_all[i]
          logacceptance = lognum - logdenom
          logu = log(runif(1))
          if (logu <= logacceptance){
            thetas[,i] = theta_new
            logtargetdensities[i] = logtarget_new
            if (!is.null(byproducts)) {byproducts[[i]] = byproduct}
          }
          else {
            # do nothing
          }
        }
        thetas_larger_sample[,((imove-1)*Ntheta+i)] = thetas[,i]
        if (!is.null(byproducts)) {byproducts_larger_sample[[((imove-1)*Ntheta+i)]] = byproducts[[i]]}
      }
      if (algorithmic_parameters$verbose) {
        Ntheta_current = Ntheta_current + Ntheta
        cat("\r Additional particles generated:",Ntheta_current,"out of",Ntheta_total,"(",floor(Ntheta_current/Ntheta_total*100),"%)")
        flush.console()
      }
    }
  }
  if (algorithmic_parameters$verbose) {
    cat("\n Done. \n")
  }
  return (list(thetas = thetas_larger_sample, byproducts = byproducts_larger_sample))
}
