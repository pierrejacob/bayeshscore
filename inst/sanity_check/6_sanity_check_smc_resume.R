##################################################################################################
# This checks that the results are saved properly and can be used to resume a run
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(19)
#--------------------------------------------------------------------------------------------
# Get model
model = get_model_lineargaussian()
#--------------------------------------------------------------------------------------------
# create some data
nobservations = 45
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
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$store_byproducts = TRUE
algorithmic_parameters$ess_threshold = 0.5
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
visualize_save = function(partial_smc_results) {
  partial_nobs = partial_smc_results$t
  thetas_smc2 = partial_smc_results$thetas_history[[partial_nobs+1]]
  normw_smc2 = partial_smc_results$normw_history[[partial_nobs+1]]
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
  plot_ESS = ggplot() + geom_line(aes(1:partial_nobs, partial_smc_results$ESS[1:partial_nobs])) + xlim(0,nobservations) + theme(legend.position="none")
  #--------------------------------------------------------------------------------------------
  # Check the log-evidence (RESCALED BY 1/t)
  results = data.frame(from = factor(rep(c("smc2"),each = partial_nobs)))
  results$time = rep(1:partial_nobs, 1)
  results$logevidence = c(cumsum(partial_smc_results$incr_logevidence[1:partial_nobs]))
  plot_logev = ggplot(results) + geom_line(aes(time, logevidence/time, color = from), size = 1) + xlim(0,nobservations) + theme(legend.position="none")
  #--------------------------------------------------------------------------------------------
  # Check the h-score (RESCALED BY 1/t)
  results$hscore = c(cumsum(partial_smc_results$incr_hscore[1:partial_nobs]))
  plot_hscore = ggplot(results) + geom_line(aes(time, hscore/time, color = from), size = 1) + xlim(0,nobservations) + theme(legend.position="none")
  #--------------------------------------------------------------------------------------------
  grid.arrange(plot_theta1, plot_ESS, plot_theta2, plot_logev,
               plot_theta3, plot_hscore, plot_theta4 , ncol = 2, nrow = 4)
}
#########################################################################################
#########################################################################################
###                                                                             #########
###  Assimilate first batch of observations: PURPOSELY INTERRUPT COMPUTATION    #########
###                                                                             #########
#########################################################################################
#########################################################################################
# Purposely assimilate the observations as two separate batch to test smc_resume features
n_first_batch = floor(nobservations/2)
#----------------------------------------------------------------------------------------
# WARNING: the save must be an RDS file (extension .rds)
algorithmic_parameters$savefilename = "partial_smc_results.rds"
# purposely set a time budget
algorithmic_parameters$time_budget = 5
#----------------------------------------------------------------------------------------
# Run SMC with a time budget
smc_result = hscore(observations[,1:n_first_batch,drop=FALSE],model,algorithmic_parameters)
# Load and visualize partial results
partial_smc_results = readRDS(algorithmic_parameters$savefilename)
visualize_save(partial_smc_results)
#----------------------------------------------------------------------------------------
# Continue to resume / interrupt run until all the observations have been assimilated
count = 0
while (partial_smc_results$t < n_first_batch) {
  # increase time budget just to make sure the while loop progresses
  algorithmic_parameters$time_budget = 2*algorithmic_parameters$time_budget
  # resume SMC run
  smc_result = smc_resume(RDSsave = partial_smc_results, new_algorithmic_parameters = algorithmic_parameters)
  # save and visualize partial results
  partial_smc_results = readRDS(algorithmic_parameters$savefilename)
  visualize_save(partial_smc_results)
}
cat("----- First batch assimilated -----\n")
#########################################################################################
#############               Assimilate newly available observations       ###############
#########################################################################################
cat("----- Assimilating second batch -----\n")
smc_result = smc_resume(RDSsave = partial_smc_results,
                        next_observations = observations[,(n_first_batch+1):nobservations,drop=FALSE],
                        new_algorithmic_parameters = algorithmic_parameters)
partial_smc_results = readRDS(algorithmic_parameters$savefilename)
visualize_save(partial_smc_results)
while (partial_smc_results$t < nobservations) {
  # increase time budget just to make sure the while loop progresses
  algorithmic_parameters$time_budget = 2*algorithmic_parameters$time_budget
  # resume SMC run
  smc_result = smc_resume(RDSsave = partial_smc_results, new_algorithmic_parameters = algorithmic_parameters)
  # save and visualize partial results
  partial_smc_results = readRDS(algorithmic_parameters$savefilename)
  visualize_save(partial_smc_results)
}
cat("----- Second batch assimilated -----\n")
