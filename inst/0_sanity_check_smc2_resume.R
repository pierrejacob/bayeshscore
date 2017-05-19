##################################################################################################
# This checks that the results are saved properly and can be used to resume a run
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(19)
module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree
#--------------------------------------------------------------------------------------------
# Get model
model = get_model_lineargaussian()
model_nolikelihood = model
model_nolikelihood$likelihood = NULL # this forces the use of SMC2
model_nolikelihood$dpredictive = NULL # this forces the use of SMC2
#--------------------------------------------------------------------------------------------
# create some data
nobservations = 20
theta_star = c(0.8,1,1,1)
sim = simulateData(model, theta = theta_star, nobservations)
X = sim$X
Y = sim$Y
observations = matrix(Y, nrow = model$dimY)
# observations in a matrix of dimensions dimY by nobservations
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$Nx = 2^5
algorithmic_parameters$Nx_max = 2^10
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$store_X = TRUE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$min_acceptance_rate = 0.3
algorithmic_parameters$nmoves = 1
algorithmic_parameters$save = TRUE
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R
#########################################################################################
#########################################################################################
#########################################################################################
#############                 DIAGNOSTICS of partial save                 ###############
#########################################################################################
#########################################################################################
#########################################################################################
########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################

# define some temporary diagnostic function
visualize_save = function(partial_smc2_results) {
  partial_nobs = partial_smc2_results$t
  thetas_smc2 = partial_smc2_results$thetas_history[[partial_nobs+1]]
  normw_smc2 = partial_smc2_results$normw_history[[partial_nobs+1]]
  #--------------------------------------------------------------------------------------------
  # Checking sample from the posterior distribution (marginal histogram)
  Ntheta = algorithmic_parameters$Ntheta
  post = data.frame(from = factor(rep(c("smc2"),each = Ntheta)))
  post$theta1 = c(thetas_smc2[1,])
  post$theta2 = c(thetas_smc2[2,])
  post$theta3 = c(thetas_smc2[3,])
  post$theta4 = c(thetas_smc2[4,])
  post$weight = c(normw_smc2)
  # plot posterior marginals
  plot_theta1 = ggplot(post) + geom_density(aes(theta1, weight = weight, fill = from), alpha = 0.6) + theme(legend.position="none")
  plot_theta2 = ggplot(post) + geom_density(aes(theta2, weight = weight, fill = from), alpha = 0.6) + theme(legend.position="none")
  plot_theta3 = ggplot(post) + geom_density(aes(theta3, weight = weight, fill = from), alpha = 0.6) + theme(legend.position="none")
  plot_theta4 = ggplot(post) + geom_density(aes(theta4, weight = weight, fill = from), alpha = 0.6) + theme(legend.position="none")
  #--------------------------------------------------------------------------------------------
  # Check ESS
  plot_ESS = ggplot() + geom_line(aes(1:partial_nobs, partial_smc2_results$ESS[1:partial_nobs])) + xlim(0,nobservations) + theme(legend.position="none")
  #--------------------------------------------------------------------------------------------
  # Check the log-evidence (RESCALED BY 1/t)
  results = data.frame(from = factor(rep(c("smc2"),each = partial_nobs)))
  results$time = rep(1:partial_nobs, 1)
  results$logevidence = c(cumsum(partial_smc2_results$incr_logevidence[1:partial_nobs]))
  plot_logev = ggplot(results) + geom_line(aes(time, logevidence/time, color = from), size = 1) + xlim(0,nobservations) + theme(legend.position="none")
  #--------------------------------------------------------------------------------------------
  # Check the h-score (RESCALED BY 1/t)
  results$hscore = c(cumsum(partial_smc2_results$incr_hscore[1:partial_nobs]))
  plot_hscore = ggplot(results) + geom_line(aes(time, hscore/time, color = from), size = 1) + xlim(0,nobservations) + theme(legend.position="none")
  #--------------------------------------------------------------------------------------------
  # estimate (marginal) filtering means from SMC2
  xmean = rep(0,partial_nobs)
  Ntheta = partial_smc2_results$algorithmic_parameters$Ntheta
  for (t in 1:partial_nobs) {
    PF = partial_smc2_results$PF_history[[t+1]]
    normw = partial_smc2_results$normw_history[[t+1]]
    xmean[t] = sum(sapply(1:Ntheta,function(i)sum(PF[[i]]$X*PF[[i]]$xnormW)*normw[i]))
  }
  # estimate (marginal) filtering means via KF (using the same particles thetas from SMC2)
  xmeanKF = rep(0,partial_nobs)
  for (t in 1:partial_nobs) {
    normw = partial_smc2_results$normw_history[[t+1]]
    thetas = partial_smc2_results$thetas_history[[t+1]]
    for (i in 1:Ntheta){
      phi = thetas[1,i]
      sigmaW2 = thetas[2,i]
      psi = thetas[3,i]
      sigmaV2 = thetas[4,i]
      initial_mean = 0
      initial_var = (sigmaW2)/(1-phi^2)
      KF = KF_filtering(partial_smc2_results$observations[,1:t,drop=FALSE],phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var)
      xmeanKF[t] = xmeanKF[t] + KF[[t]]$muX_t_t*normw[i]
    }
  }
  # plot filtering means
  plot_PF = ggplot(data.frame(time = rep(1:partial_nobs,2), marginalfiltermean = c(xmean,xmeanKF),
                         from = factor(rep(c("smc2","KF"),each=partial_nobs)))) +
    geom_line(aes(time, marginalfiltermean, color = from, linetype = from),size=1) +
    geom_point(aes(time, marginalfiltermean, color = from, shape = from),size=3) +
    xlim(0,nobservations) + theme(legend.position="none")
  grid.arrange(plot_theta1, plot_ESS, plot_theta2, plot_logev,
               plot_theta3, plot_hscore, plot_theta4, plot_PF, ncol = 2, nrow = 4)
}
#########################################################################################
#########################################################################################
#########################################################################################
#############               PURPOSELY INTERRUPT COMPUTATION               ###############
#########################################################################################
#########################################################################################
#########################################################################################
algorithmic_parameters$savefilename = "partial_smc2_results.rds"
#
# WARNING: the save must be an RDS file (extension .rds)
#
time_budget = 10
cat(strftime(Sys.time()),", alloted time =",time_budget,"sec\n")
setTimeLimit(elapsed = time_budget) # purposely set a time budget
tryCatch(hscore(observations,model_nolikelihood,algorithmic_parameters), error = function(e) cat("--- Time out ---\n"))
setTimeLimit() # This resets the time limit to infinity
partial_smc2_results = readRDS(algorithmic_parameters$savefilename)
visualize_save(partial_smc2_results)
#
# Continue to resume / interrupt run until all the observations have been assimilated
count = 1
while (partial_smc2_results$t < nobservations) {
  time_budget = 10*(2^count)
  cat(strftime(Sys.time()),", alloted time =",time_budget,"sec\n")
  setTimeLimit(elapsed = time_budget) # purposely set a time budget
  tryCatch(smc2_resume(partial_smc2_results), error = function(e) cat("--- Time out ---\n"))
  setTimeLimit() # This resets the time limit to infinity
  partial_smc2_results = readRDS(algorithmic_parameters$savefilename)
  visualize_save(partial_smc2_results)
  count = count + 1
}
setTimeLimit() # This resets the time limit to infinity

