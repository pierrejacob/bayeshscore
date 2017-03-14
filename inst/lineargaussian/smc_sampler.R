rm(list = ls())
library(HyvarinenSSM)
library(doMC)
set.seed(17)
set_global_path()
assimilate_one <- function(thetas, KFs, t, observations, model, Ntheta, ess_objective, nmoves, resampling, logtargetdensities, logw, normw){
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
      maxlogw <- max(logw_);  w <- exp(logw_ - maxlogw); normw <- w / sum(w)
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
        gamma <- seach_gamma(current_gamma, ess_given_gamma, objective = ess_objective)$x
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
    cat("Step", t, ", gamma = ", gamma, ", ESS = ", 1/(sum(normw^2)), "\n")
    if (gamma<1){
      # we need to resample and move
      # resampling step
      covariance = cov.wt(thetas, wt = normw, method = "ML")
      mean_t = covariance$center
      cov_t = matrix(covariance$cov,nrow = model$dimtheta) + diag(rep(10^(-4)/model$dimtheta),model$dimtheta)
      #(increased a little bit the diagonal to prevent degeneracy effects)
      resampled_index = resampling(normw)
      thetas = thetas[resampled_index,,drop=FALSE]
      logtargetdensities <- logtargetdensities[resampled_index]
      logw_incremental <- logw_incremental[resampled_index]
      logw <- rep(0, Ntheta)
      normw <- rep(1/Ntheta, Ntheta)
      #
      if (nmoves > 0){
        for (imove in 1:nmoves){
          theta_new_all = fast_rmvnorm(Ntheta,mean_t,cov_t)
          proposal_density_new_all <- fast_dmvnorm(theta_new_all, mean_t, cov_t)
          proposal_density_current <- fast_dmvnorm(thetas, mean_t, cov_t)
          accepts <- 0
          for (i in 1:Ntheta) {
            theta_new <- theta_new_all[i,]
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
              lognum <- logprior_theta_new + loglikelihood_new + proposal_density_current[i]
              logdenom <- logtargetdensities[i] + proposal_density_new_all[i]
              logacceptance <- lognum - logdenom
              logu <- log(runif(1))
              if (logu <= logacceptance){
                accepts <- accepts + 1
                thetas[i,] <- theta_new
                logtargetdensities[i] <- logprior_theta_new + loglikelihood_new
                logw_incremental[i] <- logw_incremental_new
                KFs[[i]] = Kalman
              }
              else {
                # do nothing
              }
            }
          }
          cat("Acceptance rate (independent proposal): ", 100*accepts/Ntheta, "%\n")
        }
      }
    }
  }
  return(list(KFs = KFs, thetas = thetas, normw = normw, logw = logw, logtargetdensities = logtargetdensities, logcst = logcst))
}


# load data
nobservations <- 50
model <- get_model_simplerlineargaussian()
theta_star <- model$theta
sim = simulateData(model, theta = theta_star, nobservations)
X = sim$X
Y = sim$Y
# observations in a matrix of dimensions dimy x nobservations
observations <- matrix(Y, nrow = model$dimY)


kalman_module <<- Module( "kalman_mod", PACKAGE = "HyvarinenSSM")

smc_sampler <- function(observations, model, algorithmic_parameters){
  Ntheta <- algorithmic_parameters$Ntheta
  nmoves = algorithmic_parameters$nmoves
  resampling = algorithmic_parameters$resampling

  ess_objective <- algorithmic_parameters$ess_threshold*algorithmic_parameters$Ntheta
  nobservations <- ncol(observations)
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
  # initialize Kalman filters
  KFs <- list()

  for (itheta in 1:Ntheta){
    theta <- thetas[itheta,]
    Kalman <- new(kalman_module$Kalman)
    Kalman$set_parameters(list(rho = theta[1], sigma = sqrt(theta[2]), eta = 1, tau = sqrt(0.8)))
    Kalman$set_observations(matrix(observations, ncol = 1))
    Kalman$first_step()
    KFs[[itheta]] <- Kalman
  }
  for (t in 1:nobservations){
    results <- assimilate_one(thetas, KFs, t, observations, model, Ntheta, ess_objective, nmoves, resampling, logtargetdensities, logw, normw)
    thetas <- results$thetas
    normw <- results$normw
    logw <- results$logw
    KFs <- results$KFs
    logtargetdensities <- results$logtargetdensities
    logevidence[t] <- results$logcst
    thetas_history[[t+1]] <- thetas
    normw_history[[t+1]] <- normw
  }
  return(list(thetas_history = thetas_history, normw_history = normw_history, logevidence = logevidence, logtargetdensities = logtargetdensities))
}

algorithmic_parameters <- list(Ntheta = 1024, resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
                               ess_threshold = 0.5, nmoves = 10)

smc_results <- smc_sampler(observations, model, algorithmic_parameters)
thetas <- smc_results$thetas_history[[nobservations+1]]
normw <- smc_results$normw_history[[nobservations+1]]

qplot(x = thetas[,1], weight = normw, geom = "blank") + geom_histogram(aes(y = ..density..))
qplot(x = thetas[,1], y = thetas[,2], alpha = normw, geom = "point") + geom_point(aes(x = theta_star[1], y = theta_star[2]), colour = "red", size = 5)



# qplot(x = thetas[,3], y = thetas[,4], alpha = normw, geom = "point")  + geom_point(aes(x = theta_star[3], y = theta_star[4]), colour = "red", size = 5)

# qplot(x = thetas[,2], weight = normw, geom = "blank") + geom_histogram(aes(y = ..density..))
# qplot(x = thetas[,3], weight = normw, geom = "blank") + geom_histogram(aes(y = ..density..))
# qplot(x = thetas[,4], weight = normw, geom = "blank") + geom_histogram(aes(y = ..density..))

# theta <- theta_star

# # calculate log-likelihood for various parameters (changing only the first component)
# library(dlm)
# phi <- theta[1]
# psi = theta[2]
# sigmaV2 = theta[3]
# sigmaW2 <- theta[4]
#
# mod <- dlm(list(m0 = 0, C0 = sigmaW2/(1-phi^2), FF = psi, V = sigmaV2, GG = phi, W = sigmaW2))
# -dlmLL(observations[1,], mod)
# library(astsa)
# astsa_results <- Kfilter0(nobservations, observations[1,], psi, 0, sigmaW2/(1-phi^2), phi, sqrt(sigmaW2), sqrt(sigmaV2))
# -astsa_results$like
#
# kalman_module <- Module( "kalman_mod", PACKAGE = "HyvarinenSSM")
#
# kalman_ll <- function(theta){
#   phi <- theta[1]
#   psi = theta[2]
#   sigmaV2 = theta[3]
#   sigmaW2 <- theta[4]
#   LGModel <- new(kalman_module$LinearGaussian)
#   LGModel$setLinearGaussianMatrices()
#   LGModel$set_parameters(list(rho = phi, sigma = sqrt(sigmaW2), eta = psi, tau = sqrt(sigmaV2)))
#   LGModel$set_observations(matrix(observations, ncol = 1))
#   Kalman <- new(kalman_module$Kalman)
#   Kalman$setLinearGaussian(LGModel)
#   Kalman$filtering()
#   return(Kalman$get_incremental_ll())
# }
#
# kalman_ll(theta)
#
#
# nparam <- 10
# theta1grid <- seq(from = 0.3, to = 0.9, length.out = nparam)
# log_p_y_hat <- c()
# log_p_y <- c()
# for (i in 1:nparam){
#   theta[1] <- theta1grid[i]
#   #
#   log_p_y <- c(log_p_y, sum(kalman_ll(theta)))
#   #
#   pf_results <- bootstrap_particle_filter(observations, model, theta, algorithmc_parameters)
#   log_p_y_hat <- c(log_p_y_hat, pf_results$log_p_y_hat)
# }
# plot(theta1grid, log_p_y_hat)
# lines(theta1grid, log_p_y)
#
# kalman_ll(theta_star)
#
