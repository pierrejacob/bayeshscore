rm(list = ls())
library(HyvarinenSSM)
library(doMC)
set.seed(19)
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
  return(list(thetas_history = thetas_history, normw_history = normw_history,
              logevidence = cumsum(logevidence), logtargetdensities = logtargetdensities))
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

algorithmic_parameters <- list(Ntheta = 1024, resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
                               ess_threshold = 0.5, nmoves = 2)
algorithmic_parameters$Nx = NULL # this uses adaptive Nx, starting with Nx = 128
algorithmic_parameters$min_acceptance_rate = 0.45


kalman_module <<- Module( "kalman_mod", PACKAGE = "HyvarinenSSM")

smc_results <- smc_sampler(observations, model, algorithmic_parameters)

###############################################################################################
###############################################################################################

# #### Sanity check: sample posterior via naive MH
# # Computes the posterior density (target)
# psi = model$psi
# sigmaV2 = model$sigmaV2
# dpost = function(theta, log = TRUE){
#   Kalman <- new(kalman_module$Kalman)
#   Kalman$set_parameters(list(rho = theta[1], sigma = sqrt(theta[2]), eta = psi, tau = sqrt(sigmaV2)))
#   Kalman$set_observations(matrix(observations, ncol = 1))
#   Kalman$first_step()
#   for (t in 1:nobservations){
#     Kalman$filtering_step(t-1)
#   }
#   loglikelihood = sum(Kalman$get_incremental_ll())
#   if (log){
#     return (model$dprior(theta,log = TRUE) + loglikelihood)
#   } else {
#     return (exp((model$dprior(theta,log = TRUE) + loglikelihood)))
#   }
# }
#
# MH_cov = (cov.wt(thetas_smc,wt=normw_smc)$cov)/5
# M = 10000
# burnin = 9000
# thetas_MH = matrix(NA,ncol = 2,nrow = M)
# ###
# print(paste("Started at:",Sys.time()))
# progbar = txtProgressBar(min = 0,max = M,style=3)
# count = 1
# time_start = proc.time()
# ###
# thetas_MH[1,] = c(0.7,1)
# accepts = 0
# ###
# for (i in 2:M){
#   theta_new = fast_rmvnorm(1,thetas_MH[i-1,],MH_cov)
#   if (model$dprior(theta_new) == -Inf){
#     thetas_MH[i,] = thetas_MH[i-1,]
#     count = count + 1
#     setTxtProgressBar(progbar, count)
#     next
#   } else {
#     logacceptance = dpost(theta_new) - dpost(thetas_MH[i-1,])
#     logu = log(runif(1))
#     if (logu <= logacceptance){
#       accepts = accepts + 1
#       thetas_MH[i,] = theta_new
#     } else {
#       thetas_MH[i,] = thetas_MH[i-1,]
#     }
#   }
#   ###
#   count = count + 1
#   setTxtProgressBar(progbar, count)
# }
# time_end = proc.time()-time_start
# print(time_end)
# cat("acceptance rate = ", accepts/(M-1))
# # Traceplots
# par(mfrow=c(2,1))
# index = seq(burnin,M,by = 1)
# plot(index,thetas_MH[index,1],type='l')
# plot(index,thetas_MH[index,2],type='l')
# par(mfrow=c(1,1))
#
#
#
# ### Check SMC_2 output
# module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
# TreeClass <<- module_tree$Tree
# algorithmic_parameters$progress = TRUE
# smc2_results <- hscore_continuous(observations, model, algorithmic_parameters)
# thetas_smc2 <- smc2_results$thetas
# normw_smc2 <- smc2_results$thetanormw
#
# # Visual (qualitative) diagnostic
# M = 100
# x = seq(0.01, 1.5, length.out = M)
# y = seq(0.5,10, length.out = M)
# z = matrix(NA,ncol = M,nrow = M)
# print(paste("Started at:",Sys.time()))
# progbar = txtProgressBar(min = 0,max = M*M,style=3)
# count = 0
# time_start = proc.time()
# for (i in 1:M){
#   for (j in 1:M){
#     z[i,j] = dpost(c(x[i],y[j]))
#     count = count + 1
#     setTxtProgressBar(progbar, count)
#   }
# }
# time_end = proc.time()-time_start
# print(time_end)
# contour_df = expand.grid(x = x, y = y)
# contour_df$z = c(z)
#
#
# # Checking sample from the posterior distribution
# ### BLACK = SMC output
# ### YELLOW = MH output
# ### PURPLE = SMC_2 output
# ggplot() +
#   geom_point(aes(thetas_smc[,1], thetas_smc[,2], alpha = normw_smc)) +
#   geom_contour(data = contour_df,aes(x,y,z=z,color = ..level..)) +
#   scale_colour_gradient(low="black", high="red") +
#   geom_point(aes(thetas_MH[index,1],thetas_MH[index,2]),color="yellow", size = 1, shape = 3) +
#   geom_point(aes(thetas_smc2[,1], thetas_smc2[,2], alpha = normw_smc2),color="purple", shape=15) +
#   geom_point(aes(x = theta_star[1], y = theta_star[2]), colour = "red", size = 10)


### Check SMC_2 output
module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree
algorithmic_parameters$progress = TRUE
algorithmic_parameters$store = TRUE

# Tempered SMC^2
smc2_results_temp <- hscore_continuous(observations, model, algorithmic_parameters)
# Regular SMC^2
smc2_results <- hscore_continuous_no_tempering(observations, model, algorithmic_parameters)


########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################
time_t = 50
#
thetas_smc <- smc_results$thetas_history[[time_t+1]]
normw_smc <- smc_results$normw_history[[time_t+1]]
#
thetas_smc2 <- smc2_results$thetas_history[[time_t+1]]
normw_smc2 <- smc2_results$weights_history[[time_t+1]]
#
thetas_smc2_temp <- smc2_results_temp$thetas_history[[time_t+1]]
normw_smc2_temp <- smc2_results_temp$normw_history[[time_t+1]]

# Checking sample from the posterior distribution
library(gridExtra)
plot_theta1 = ggplot() +
  geom_density(aes(thetas_smc[,1], weight = normw_smc, alpha = 0.1), fill = "green") +
  geom_density(aes(thetas_smc2[,1], weight = normw_smc2, alpha = 0.1), fill = "purple") +
  geom_density(aes(thetas_smc2_temp[,1], weight = normw_smc2_temp, alpha = 0.1), fill = "blue")
# geom_density(aes(thetas_MH[,1], alpha = 0.1), fill = "yellow")
plot_theta2 = ggplot() +
  geom_density(aes(thetas_smc[,2], weight = normw_smc, alpha = 0.1), fill = "green") +
  geom_density(aes(thetas_smc2[,2], weight = normw_smc2, alpha = 0.1), fill = "purple") +
  geom_density(aes(thetas_smc2_temp[,2], weight = normw_smc2_temp, alpha = 0.1), fill = "blue")
# geom_density(aes(thetas_MH[,2], alpha = 0.1), fill = "yellow")
grid.arrange(plot_theta1, plot_theta2, ncol = 2)

# Check the log-evidence (RESCALED BY 1/t)
ggplot() +
  geom_line(aes(1:time_t,smc_results$logevidence[1:time_t]/1:time_t), color = "green", size = 1) +
  geom_line(aes(1:time_t,smc2_results$logevidence[1:time_t]/1:time_t), color = "purple", size = 1) +
  geom_line(aes(1:time_t,smc2_results_temp$logevidence[1:time_t]/1:time_t), color = "blue", size = 1)

# Check the h-score (RESCALED BY 1/t)
ggplot() +
  geom_line(aes(1:time_t,smc2_results$hscore[1:time_t]/1:time_t), color = "purple", size = 1) +
  geom_line(aes(1:time_t,smc2_results_temp$Hscore[1:time_t]/1:time_t), color = "blue", size = 1)
