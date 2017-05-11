# This function does one step of SMC by assimilating one observation
assimilate_one_smc <- function(thetas, t, observations, model, Ntheta, ess_objective, nmoves,
                               resampling, logtargetdensities, logw, normw, byproducts, verbose = FALSE){
  current_gamma <- 0
  logcst <- 0
  logw_incremental <- rep(NA, Ntheta)
  # compute the log-incremental weight using byproduct (e.g. Kalman filter) or analytically via the model
  if (!is.null(byproducts)){
    for (itheta in 1:Ntheta){
      byproduct = model$update_byproduct(byproducts[[itheta]], t, thetas[,itheta], observations)
      logw_incremental[itheta] = model$dpredictive(observations,t,thetas[,itheta],byproduct,log = TRUE)
    }
  } else {
    for (itheta in 1:Ntheta){
      logw_incremental[itheta] = model$dpredictive(observations,t,thetas[,itheta],log = TRUE)
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
    # display diagnostic
    if (verbose){
      cat("Step", t, ", gamma = ", gamma, ", ESS = ", 1/(sum(normw^2)), "\n")
    }
    if (gamma<1){
      # we need to resample and move
      # resampling step
      covariance = cov.wt(t(thetas), wt = normw, method = "ML")
      mean_t = covariance$center
      cov_t = covariance$cov + diag(rep(10^(-4)/model$dimtheta),model$dimtheta)
      # (increased a little bit the diagonal to prevent degeneracy effects)
      resampled_index = resampling(normw)
      thetas = thetas[,resampled_index,drop=FALSE]
      logtargetdensities <- logtargetdensities[resampled_index]
      logw_incremental <- logw_incremental[resampled_index]
      logw <- rep(0, Ntheta)
      normw <- rep(1/Ntheta, Ntheta)
      #
      if (nmoves > 0){
        for (imove in 1:nmoves){
          theta_new_all = fast_rmvnorm_transpose(Ntheta,mean_t,cov_t)
          log_proposal_density_new_all <- fast_dmvnorm_transpose(theta_new_all, mean_t, cov_t)
          log_proposal_density_current <- fast_dmvnorm_transpose(thetas, mean_t, cov_t)
          accepts <- 0
          for (i in 1:Ntheta) {
            theta_new <- theta_new_all[,i]
            logprior_theta_new <- model$dprior(theta_new, log = TRUE)
            if (is.infinite(logprior_theta_new)){
              next
            } else {
              # compute the log-likelihood of proposed theta using byproduct (e.g. Kalman filter) or analytically via the model
              if (!is.null(byproducts)){
                byproduct = model$initialize_byproducts(theta_new, observations, Ntheta)
                for (time_index in 1:t){
                  byproduct = model$update_byproduct(byproduct, time_index, theta_new, observations)
                }
                incremental_ll_new = model$dpredictive(observations,t,theta_new,byproduct,log = TRUE)
              } else {
                incremental_ll_new = model$dpredictive(observations,t,theta_new,log = TRUE)
              }
              loglikelihood_new <- gamma * incremental_ll_new
              logw_incremental_new <- incremental_ll_new
              if (t > 1){
                # compute the log-likelihood of proposed theta using byproduct (e.g. Kalman filter) or analytically via the model
                if (!is.null(byproducts)){
                  loglikelihood_new <- loglikelihood_new + model$likelihood(observations,t-1,theta_new,byproduct,log = TRUE)
                } else {
                  loglikelihood_new <- loglikelihood_new + model$likelihood(observations,t-1,theta_new,log = TRUE)
                }
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
                byproducts[[i]] = byproduct
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
  return(list(thetas = thetas, normw = normw, logw = logw, logtargetdensities = logtargetdensities,
              logcst = logcst, byproducts = byproducts))
}
