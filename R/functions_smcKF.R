# This function does one step of SMC by assimilating one observation
assimilate_one_smcKF <- function(thetas, KFs, t, observations, model, Ntheta, ess_objective, nmoves,
                               resampling, logtargetdensities, logw, normw, verbose = FALSE){
  current_gamma <- 0
  logcst <- 0
  logw_incremental <- rep(0, Ntheta)
  for (itheta in 1:Ntheta){
    KF <- KFs[[itheta]]
    KF$filtering_step(t-1)
    logw_incremental[itheta] <- KF$get_incremental_ll()[t]
    KFs[[itheta]] <- KF
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
              # KF on the proposed theta
              Kalman <- new(kalman_module$Kalman)
              Kalman$set_parameters(list(rho = theta_new[1], sigma = sqrt(theta_new[2]), eta = 1, tau = sqrt(0.8)))
              Kalman$set_observations(matrix(observations, ncol = 1))
              Kalman$first_step()
              for (time_index in 1:t){
                Kalman$filtering_step(time_index-1)
              }
              incremental_ll_new <- Kalman$get_incremental_ll()
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
                KFs[[i]] = Kalman
              }
              else {
                # do nothing
              }
            }
          }
          if (verbose){
            cat("Acceptance rate (independent proposal): ", 100*accepts/Ntheta, "%\n")
          }
        }
      }
    }
  }
  return(list(KFs = KFs, thetas = thetas, normw = normw, logw = logw,
              logtargetdensities = logtargetdensities, logcst = logcst))
}
