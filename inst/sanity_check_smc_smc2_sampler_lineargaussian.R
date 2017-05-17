rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(19)
module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree
#--------------------------------------------------------------------------------------------
# create data
nobservations <- 15
model <- get_model_lineargaussian()
theta_star <- c(0.8,1,1,1)
#--------------------------------------------------------------------------------------------
sim = simulateData(model, theta = theta_star, nobservations)
X = sim$X
Y = sim$Y
observations <- matrix(Y, nrow = model$dimY)
# observations in a matrix of dimensions dimY by nobservations
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters <- list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$min_acceptance_rate = 0.2
algorithmic_parameters$nmoves = 2
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R
#--------------------------------------------------------------------------------------------
model_withoutlikelihood = model
model_withoutlikelihood$likelihood = NULL # this forces the use of SMC2
model_withoutlikelihood$dpredictive = NULL # this forces the use of SMC2
#--------------------------------------------------------------------------------------------
### Run SMC_sampler
smc_results = smc_sampler(observations, model, algorithmic_parameters)
### Run SMC2_sampler
smc2_results = smc2_sampler(observations, model_withoutlikelihood, algorithmic_parameters)
### Run hscore via SMC
hsmc_results = hscore(observations, model, algorithmic_parameters)
### Run hscore via SMC2
hsmc2_results = hscore(observations, model_withoutlikelihood, algorithmic_parameters)
#--------------------------------------------------------------------------------------------
#########################################################################################
########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################
#########################################################################################
thetas_smc <- smc_results$thetas_history[[nobservations+1]]
normw_smc <- smc_results$normw_history[[nobservations+1]]
#
thetas_smc2 <- smc2_results$thetas_history[[nobservations+1]]
normw_smc2 <- smc2_results$normw_history[[nobservations+1]]
#
thetas_hsmc <- hsmc_results$thetas_history[[nobservations+1]]
normw_hsmc <- hsmc_results$normw_history[[nobservations+1]]
#
thetas_hsmc2 <- hsmc2_results$thetas_history[[nobservations+1]]
normw_hsmc2 <- hsmc2_results$normw_history[[nobservations+1]]
#--------------------------------------------------------------------------------------------
# Checking sample from the posterior distribution (marginal histogram)
Ntheta = algorithmic_parameters$Ntheta
post = data.frame(from = factor(rep(c("smc","smc2","hsmc","hsmc2"),each = Ntheta)))
post$theta1 = c(thetas_smc[1,],thetas_smc2[1,],thetas_hsmc[1,],thetas_hsmc2[1,])
post$theta2 = c(thetas_smc[2,],thetas_smc2[2,],thetas_hsmc[2,],thetas_hsmc2[2,])
post$theta3 = c(thetas_smc[3,],thetas_smc2[3,],thetas_hsmc[3,],thetas_hsmc2[3,])
post$theta4 = c(thetas_smc[4,],thetas_smc2[4,],thetas_hsmc[4,],thetas_hsmc2[4,])
post$weight = c(normw_smc,normw_smc2,normw_hsmc,normw_hsmc2)
# plot posterior marginals
plot_theta1 = ggplot(post) +  geom_density(aes(theta1, weight = weight, fill = from), alpha = 0.6) + theme(legend.position="none")
plot_theta2 = ggplot(post) +  geom_density(aes(theta2, weight = weight, fill = from), alpha = 0.6) + theme(legend.position="none")
plot_theta3 = ggplot(post) +  geom_density(aes(theta3, weight = weight, fill = from), alpha = 0.6) + theme(legend.position="none")
plot_theta4 = ggplot(post) +  geom_density(aes(theta4, weight = weight, fill = from), alpha = 0.6)
grid.arrange(plot_theta1, plot_theta2, plot_theta3, plot_theta4, ncol = 4, widths=c(1,1,1,1.5))
#--------------------------------------------------------------------------------------------
# Check the log-evidence (RESCALED BY 1/t)
results = data.frame(from = factor(rep(c("smc","smc2","hsmc","hsmc2"),each = nobservations)))
results$time = rep(1:nobservations, 4)
results$logevidence = c(smc_results$logevidence,smc2_results$logevidence,hsmc_results$logevidence,hsmc2_results$logevidence)
ggplot(results) + geom_line(aes(time, logevidence/time, color = from), size = 1)
#--------------------------------------------------------------------------------------------
