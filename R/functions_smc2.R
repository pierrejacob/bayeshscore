#######################################################################################################
# This function performs one step of particle filter for a given theta
#######################################################################################################
filter_step <- function(observationt, t, model, theta, X, xnormW, tree, resampling){
  Xnew = X
  xnormWnew = xnormW
  log_z_incremental = NA #likelihood estimate
  ancestors <- resampling(xnormW) #sample the ancestors' indexes
  X_current <- X[ancestors,]
  if (is.null(dim(X_current))) X_current <- matrix(X_current, ncol = model$dimX)
  X_current <- model$rtransition(X_current, t, theta)
  logW <- model$dobs(observationt, X_current, t, theta,log = TRUE)
  maxlogW <- max(logW)
  W <- exp(logW - maxlogW)
  log_z_incremental <- log(mean(W)) + maxlogW #udpate likelihood estimate
  xnormWnew <- W / sum(W)
  Xnew <- X_current
  tree$update(t(Xnew), ancestors-1)
  return(list(X = Xnew, xnormW = xnormWnew, log_z_incremental = log_z_incremental, tree = tree))
}

#######################################################################################################
# This function increases Nx the number of particles in x
#######################################################################################################
increase_Nx <- function(observations, t, model, thetas, PFs, Ntheta){
  Nx = PFs[[1]]$Nx
  Nx_new <- 2*Nx #by default we multiply the number of particles by 2
  # Initialize empty arrays of larger size
  X_new = array(NA,dim = c(Nx_new, model$dimX, Ntheta)) #Nx particles (most recent) for each theta (size = Nx,dimX,Ntheta)
  xnormW_new = matrix(NA, nrow = Nx_new, ncol = Ntheta) #matrix of corresponding normalized X-weights (size = Nx,Ntheta)
  log_z_new = rep(0, Ntheta) #matrix of log-likelihood estimates (size = Ntheta)
  # Construct list of trees to store paths (one tree for each theta)
  trees_new = list()
  for (i in 1:Ntheta){
    trees_new[[i]] = new(TreeClass, Nx_new, 10*Nx_new, model$dimX)
  }
  # Perform a conditional particle filter for each theta
  for (i in 1:Ntheta){
    PF = PFs[[i]]
    xnormW = PF$xnormW
    current_tree = PF$tree
    current_path = current_tree$get_path(sample(x = 0:(Nx-1), size = 1, replace = TRUE, prob = xnormW))
    CPF = conditional_particle_filter(matrix(observations[,1:t],ncol = t), model, thetas[,i], Nx_new, path = current_path)
    PFs[[i]] = CPF
  }
  return (PFs)
}

#######################################################################################################
# This function does one step of SMC^2 by assimilating the next observation
#######################################################################################################
assimilate_one_smc2 <- function(thetas, PFs, t, observations, model, Ntheta, ess_objective, nmoves,
                                resampling, logtargetdensities, logw, normw, verbose = FALSE,
                                adaptNx = FALSE, min_acceptance_rate = 0.10){
  Nx = PFs[[1]]$Nx
  current_gamma <- 0
  logcst <- 0
  trees = list()
  logw_incremental = rep(NA, Ntheta)
  rejuvenation_time = NA
  rejuvenation_accept_rate = NA
  increase_Nx_times = NA
  increase_Nx_values = NA
  if (t == 1){
    for (itheta in 1:Ntheta){
      logw_incremental[itheta] <- PFs[[itheta]]$incremental_ll[1]
      trees[[itheta]] <- PFs[[itheta]]$tree
      ######################################################################################
      # the argument PFs was initialized externally using the first observation, so it
      # already contains all the correct values
      ######################################################################################
    }
  }
  else {
    for (itheta in 1:Ntheta){
      PF = PFs[[itheta]]
      PF_next = filter_step(observations[,t], t, model, thetas[,itheta], PF$X, PF$xnormW, PF$tree, resampling)
      # perform one step of filtering and update the current particle filter for particle itheta
      logw_incremental[itheta] <- PF_next$log_z_incremental
      PF$X = PF_next$X
      PF$xnormW = PF_next$xnormW
      PF$tree = PF_next$tree
      PFs[[itheta]] <- PF
    }
  }
  while (current_gamma < 1){
    ess_given_gamma <- function(gamma){
      logw_ <- logw + (gamma - current_gamma) * logw_incremental
      maxlogw <- max(logw_)
      w <- exp(logw_ - maxlogw)
      normw <- w / sum(w)
      return(1/(sum(normw^2)))
    }
    # try gamma = 1 first
    if (ess_given_gamma(1) > ess_objective){
      gamma <- 1
    } else {
      if (ess_given_gamma(current_gamma) < ess_objective){
        gamma <- current_gamma
        print("warning! ESS at current gamma too low; something went wrong.")
      } else {
        gamma <- search_gamma(current_gamma, ess_given_gamma, objective = ess_objective)$x
      }
    }
    # now we've found our gamma
    logw_incremental_gamma <- (gamma - current_gamma) * logw_incremental
    logtargetdensities <- logtargetdensities + logw_incremental_gamma
    current_gamma <- gamma
    # compute increment to the normalizing constant
    maxlogw <- max(logw_incremental_gamma)
    w <- exp(logw_incremental_gamma - maxlogw)
    logcst <- logcst + log(sum(normw * w)) + maxlogw
    # normalize weights
    logw <- logw + logw_incremental_gamma
    w <- exp(logw - max(logw))
    normw <- w / sum(w)
    ##
    if (verbose){
      cat("Step", t, ", gamma = ", gamma, ", ESS = ", 1/(sum(normw^2)), "\n")
    }
    if (gamma<1){
      # we need to resample and move
      # resampling step
      covariance = cov.wt(t(thetas), wt = normw, method = "ML")
      mean_t = covariance$center
      cov_t = matrix(covariance$cov,nrow = model$dimtheta) + diag(rep(10^(-4)/model$dimtheta),model$dimtheta)
      #(increased a little bit the diagonal to prevent degeneracy effects)
      resampled_index = resampling(normw)
      thetas = thetas[,resampled_index,drop=FALSE]
      logtargetdensities <- logtargetdensities[resampled_index]
      logw_incremental <- logw_incremental[resampled_index]
      logw <- rep(0, Ntheta)
      normw <- rep(1/Ntheta, Ntheta)
      #
      if (nmoves > 0){
        rejuvenation_time = t
        for (imove in 1:nmoves){
          theta_new_all = fast_rmvnorm(Ntheta,mean_t,cov_t)
          log_proposal_density_new_all <- fast_dmvnorm(theta_new_all, mean_t, cov_t)
          log_proposal_density_current <- fast_dmvnorm(thetas, mean_t, cov_t)
          accepts <- 0
          for (i in 1:Ntheta) {
            theta_new <- theta_new_all[,i]
            logprior_theta_new <- model$dprior(theta_new, log = TRUE)
            if (is.infinite(logprior_theta_new)){
              next
            } else {
              #the CPF function performs a regular PF when no conditioning path is provided
              PF <- conditional_particle_filter(matrix(observations[,1:t],ncol = t), model, theta_new, Nx)
              incremental_ll_new <- PF$incremental_ll
              loglikelihood_new <- gamma * incremental_ll_new[t]
              logw_incremental_new <- incremental_ll_new[t]
              if (t > 1){
                loglikelihood_new <- loglikelihood_new + sum(incremental_ll_new[1:(t-1)])
              }
              logtarget_new = logprior_theta_new + loglikelihood_new
              lognum <- logtarget_new + log_proposal_density_current[i]
              logdenom <- logtargetdensities[i] + log_proposal_density_new_all[i]
              logacceptance <- lognum - logdenom
              logu <- log(runif(1))
              if (logu <= logacceptance){
                accepts <- accepts + 1
                thetas[,i] <- theta_new
                logtargetdensities[i] <- logtarget_new
                logw_incremental[i] <- logw_incremental_new
                PFs[[i]] = PF
              }
              else {
                # do nothing
              }
            }
          }
          rejuvenation_accept_rate = accepts/Ntheta
          if (verbose){
            cat("Acceptance rate (independent proposal): ", 100*rejuvenation_accept_rate, "%\n")
          }
          if (adaptNx){
            if (rejuvenation_accept_rate < min_acceptance_rate){
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
  return(list(PFs = PFs, thetas = thetas, normw = normw, logw = logw,
              logtargetdensities = logtargetdensities, logcst = logcst,
              rejuvenation_time = rejuvenation_time, rejuvenation_accept_rate = rejuvenation_accept_rate,
              increase_Nx_times = increase_Nx_times, increase_Nx_values = increase_Nx_values))
}

