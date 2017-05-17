#######################################################################################################
# This function performs one step of particle filter for a given theta
#######################################################################################################
filter_step = function(observationt, t, model, theta, X, xnormW, tree, resampling){
  log_z_incremental = NA #likelihood estimate
  ancestors = resampling(xnormW) #sample the ancestors' indexes
  X_current = X[,ancestors,drop=FALSE]
  Xnew = model$rtransition(X_current, t, theta)
  logW = model$dobs(observationt, Xnew, t, theta,log = TRUE)
  maxlogW = max(logW)
  W = exp(logW - maxlogW)
  log_z_incremental = log(mean(W)) + maxlogW #udpate likelihood estimate
  xnormWnew = W / sum(W)
  tree$update(Xnew, ancestors-1)
  return(list(X = Xnew, xnormW = xnormWnew, log_z_incremental = log_z_incremental, tree = tree))
}

#######################################################################################################
# This function increases Nx the number of particles in x
#######################################################################################################
increase_Nx = function(observations, t, model, thetas, PFs, Ntheta){
  Nx = PFs[[1]]$Nx
  Nx_new = 2*Nx #by default we multiply the number of particles by 2
  # Perform a conditional particle filter for each theta
  for (i in 1:Ntheta){
    current_path = (PFs[[i]]$tree)$get_path(sample(x = 0:(Nx-1), size = 1, replace = TRUE, prob = PFs[[i]]$xnormW))
    CPF = conditional_particle_filter(observations[,1:t,drop=FALSE], model, thetas[,i], Nx_new, current_path)
    PFs[[i]] = CPF
  }
  return (PFs)
}

#######################################################################################################
# This function does one step of SMC2 by assimilating the next observation
#######################################################################################################
assimilate_one_smc2 = function(thetas, PFs, t, observations, model,
                               logtargetdensities, logw, normw, algorithmic_parameters){
  # Parse algorithmic parameters and set flags accordingly
  Ntheta = algorithmic_parameters$Ntheta
  nmoves = algorithmic_parameters$nmoves
  resampling = algorithmic_parameters$resampling
  adaptNx = algorithmic_parameters$adaptNx
  min_acceptance_rate = algorithmic_parameters$min_acceptance_rate
  ess_objective = algorithmic_parameters$ess_threshold*algorithmic_parameters$Ntheta
  verbose = algorithmic_parameters$verbose
  # initialize variables
  Nx = PFs[[1]]$Nx
  current_gamma = 0
  logcst = 0
  logw_incremental = rep(NA, Ntheta)
  rejuvenation_time = NA
  rejuvenation_rate = NA
  increase_Nx_times = NA
  increase_Nx_values = NA
  if (t == 1){
    for (i in 1:Ntheta){
      logw_incremental[i] = PFs[[i]]$incremental_ll[1]
      ######################################################################################
      # the argument PFs was initialized externally using the first observation, so it
      # already contains all the correct values
      ######################################################################################
    }
  }
  else {
    for (i in 1:Ntheta){
      PF_next = filter_step(observations[,t,drop=FALSE],t,model,thetas[,i],PFs[[i]]$X,PFs[[i]]$xnormW,PFs[[i]]$tree,resampling)
      # perform one step of filtering and update the current particle filter for particle theta i
      logw_incremental[i] = PF_next$log_z_incremental
      PFs[[i]]$X = PF_next$X
      PFs[[i]]$xnormW = PF_next$xnormW
      PFs[[i]]$tree = PF_next$tree
    }
  }
  while (current_gamma < 1){
    ess_given_gamma = function(gamma){
      logw_ = logw + (gamma - current_gamma) * logw_incremental
      maxlogw = max(logw_)
      w = exp(logw_ - maxlogw)
      normw = w / sum(w)
      return(1/(sum(normw^2)))
    }
    # try gamma = 1 first
    if (ess_given_gamma(1) > ess_objective){
      gamma = 1
    } else {
      if (ess_given_gamma(current_gamma) < ess_objective){
        gamma = current_gamma
        cat(">>>>>> WARNING: ESS at current gamma too low; something went wrong <<<<<<\n")
      } else {
        gamma = search_gamma(current_gamma, ess_given_gamma, objective = ess_objective)$x
      }
    }
    # now we've found our gamma
    logw_incremental_gamma = (gamma - current_gamma) * logw_incremental
    logtargetdensities = logtargetdensities + logw_incremental_gamma
    current_gamma = gamma
    # compute increment to the normalizing constant
    maxlogw = max(logw_incremental_gamma)
    w = exp(logw_incremental_gamma - maxlogw)
    logcst = logcst + log(sum(normw * w)) + maxlogw
    # normalize weights
    logw = logw + logw_incremental_gamma
    w = exp(logw - max(logw))
    normw = w / sum(w)
    ESS = 1/(sum(normw^2))
    # display diagnostic
    if (verbose){
      cat("Step", t, ", gamma = ", gamma, ", ESS = ", ESS, "\n")
    }
    if (gamma<1){
      # we need to resample and move
      # First we get the proposal for the move steps. Note: proposalmove is a list with the following fields:
      # proposalmove$r : sampler from the proposal
      # proposalmove$d : corresponding density function
      # see set_default_algorithmic_parameters (util_default.R) for more details
      proposalmove = algorithmic_parameters$proposalmove(thetas,normw,model)
      # Resample particles theta according to normalized weights
      resampled_index = resampling(normw)
      thetas = thetas[,resampled_index,drop=FALSE]
      PFs = PFs[resampled_index]
      logtargetdensities = logtargetdensities[resampled_index]
      logw_incremental = logw_incremental[resampled_index]
      logw = rep(0, Ntheta)
      normw = rep(1/Ntheta, Ntheta)
      #
      if (nmoves > 0){
        rejuvenation_time = t
        for (imove in 1:nmoves){
          theta_new_all = proposalmove$r(Ntheta)
          log_proposal_density_new_all = proposalmove$d(theta_new_all,log=TRUE)
          log_proposal_density_current = proposalmove$d(thetas,log=TRUE)
          accepts = 0
          for (i in 1:Ntheta) {
            theta_new = theta_new_all[,i]
            logprior_theta_new = model$dprior(theta_new, log = TRUE)
            if (is.infinite(logprior_theta_new)){
              next
            } else {
              #the CPF function performs a regular PF when no conditioning path is provided
              PF_new = conditional_particle_filter(observations[,1:t,drop=FALSE], model, theta_new, Nx)
              incremental_ll_new = PF_new$incremental_ll
              loglikelihood_new = gamma * incremental_ll_new[t]
              logw_incremental_new = incremental_ll_new[t]
              if (t > 1){
                loglikelihood_new = loglikelihood_new + sum(incremental_ll_new[1:(t-1)])
              }
              logtarget_new = logprior_theta_new + loglikelihood_new
              lognum = logtarget_new + log_proposal_density_current[i]
              logdenom = logtargetdensities[i] + log_proposal_density_new_all[i]
              logacceptance = lognum - logdenom
              logu = log(runif(1))
              if (logu <= logacceptance){
                accepts = accepts + 1
                thetas[,i] = theta_new
                logtargetdensities[i] = logtarget_new
                logw_incremental[i] = logw_incremental_new
                PFs[[i]] = PF_new
              }
              # otherwise do nothing (i.e. keep the current particles)
            }
          }
          rejuvenation_rate = accepts/Ntheta
          if (verbose){
            cat("Acceptance rate (independent proposal): ", 100*rejuvenation_rate, "%\n")
          }
          if (adaptNx){
            Nx_new = 2*(PFs[[1]]$Nx)
            if ((Nx_new <= algorithmic_parameters$Nx_max)&&(rejuvenation_rate < min_acceptance_rate)){
              # Increase the number Nx of particles for each theta
              PFs = increase_Nx(observations, t, model, thetas, PFs, Ntheta)
              Nx = PFs[[1]]$Nx
              increase_Nx_times = t
              increase_Nx_values = PFs[[1]]$Nx
              if (verbose){
                cat("Nx increased to: ", increase_Nx_values, "\n")
              }
            }
          }
        }
      }
    }
  }
  return(list(thetas = thetas, normw = normw, logw = logw, logtargetdensities = logtargetdensities,
              logcst = logcst, PFs = PFs, rejuvenation_time = rejuvenation_time,
              rejuvenation_rate = rejuvenation_rate, ESS = ESS,
              increase_Nx_times = increase_Nx_times, increase_Nx_values = increase_Nx_values))
}

