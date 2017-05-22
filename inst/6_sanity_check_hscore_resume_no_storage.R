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
nobservations = 100
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
  grid.arrange(plot_logev, plot_hscore,ncol = 1, nrow = 2)
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
previous_t = 0
while (partial_smc_results$t < n_first_batch) {
  # increase time budget just to make sure the while loop progresses
  if (previous_t==partial_smc_results$t) {algorithmic_parameters$time_budget = 2*algorithmic_parameters$time_budget}
  previous_t = partial_smc_results$t
  # resume SMC run
  smc_result = hscore_resume(RDSsave = partial_smc_results, new_algorithmic_parameters = algorithmic_parameters)
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
previous_t = 0
while (partial_smc_results$t < nobservations) {
  # increase time budget just to make sure the while loop progresses
  if (previous_t==partial_smc_results$t) {algorithmic_parameters$time_budget = 2*algorithmic_parameters$time_budget}
  previous_t = partial_smc_results$t
  # resume SMC run
  smc_result = hscore_resume(RDSsave = partial_smc_results, new_algorithmic_parameters = algorithmic_parameters)
  # save and visualize partial results
  partial_smc_results = readRDS(algorithmic_parameters$savefilename)
  visualize_save(partial_smc_results)
}
cat("----- Second batch assimilated -----\n")
