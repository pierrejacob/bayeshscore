rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(9) #the hscore goes wrong after the rejuvenation steps (e.g. Nx = 2^6, Ntheta = 2^10)
# however, things work fine when the rejuvenation triggers the increase Nx step .... weird ....
# And also: when no rejuvenation involved, things work fine ...
# e.g. CHANGE THE MINIMUM ACCEPTANCE RATE from 0.2 (fails) to 0.3 (works) TO SEE WHAT HAPPENS
# ... the hscore suddenly JUMPS every time there is a rejuvenation without Nx increase ... WHY ?????

#--------------------------------------------------------------------------------------------
# create data
nobservations <- 15
# model <- get_model_simplerlineargaussian()
# theta_star <- c(0.8,1,model$psi,model$sigmaV2)
model <- get_model_lineargaussian()
theta_star <- c(0.8,1,1,1)

sim = simulateData(model, theta = theta_star, nobservations)
X = sim$X
Y = sim$Y
observations <- matrix(Y, nrow = model$dimY)
# observations in a matrix of dimensions dimY by nobservations

#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters <- list()
algorithmic_parameters$Ntheta = 2^4
algorithmic_parameters$Nx = 2^5
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$store_X = TRUE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$min_acceptance_rate = 0.2
algorithmic_parameters$nmoves = 2
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R
#--------------------------------------------------------------------------------------------
### Run SMC
smc_results = hscore(observations, model, algorithmic_parameters)
### Run SMC_2
module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree
model_withoutlikelihood = model
model_withoutlikelihood$likelihood = NULL # this forces the use of SMC2
model_withoutlikelihood$dpredictive = NULL # this forces the use of SMC2
smc2_results = hscore(observations, model_withoutlikelihood, algorithmic_parameters)


#--------------------------------------------------------------------------------------------
#########################################################################################
########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################
thetas_smc <- smc_results$thetas_history[[nobservations+1]]
normw_smc <- smc_results$normw_history[[nobservations+1]]
#
thetas_smc2 <- smc2_results$thetas_history[[nobservations+1]]
normw_smc2 <- smc2_results$normw_history[[nobservations+1]]

#--------------------------------------------------------------------------------------------
# Checking sample from the posterior distribution (marginal histogram)
Ntheta = algorithmic_parameters$Ntheta
post = data.frame(from = factor(rep(c("smc","smc2"),each = Ntheta)))
post$theta1 = c(thetas_smc[1,],thetas_smc2[1,])
post$theta2 = c(thetas_smc[2,],thetas_smc2[2,])
post$theta3 = c(thetas_smc[3,],thetas_smc2[3,])
post$theta4 = c(thetas_smc[4,],thetas_smc2[4,])
post$weight = c(normw_smc,normw_smc2)
# plot posterior marginals
plot_theta1 = ggplot(post) +  geom_density(aes(theta1, weight = weight, fill = from), alpha = 0.6) + theme(legend.position="none")
plot_theta2 = ggplot(post) +  geom_density(aes(theta2, weight = weight, fill = from), alpha = 0.6) + theme(legend.position="none")
plot_theta3 = ggplot(post) +  geom_density(aes(theta3, weight = weight, fill = from), alpha = 0.6) + theme(legend.position="none")
plot_theta4 = ggplot(post) +  geom_density(aes(theta4, weight = weight, fill = from), alpha = 0.6)
grid.arrange(plot_theta1, plot_theta2, plot_theta3, plot_theta4, ncol = 4, widths=c(1,1,1,1.5))

#--------------------------------------------------------------------------------------------
# Check the log-evidence (RESCALED BY 1/t)
results = data.frame(from = factor(rep(c("smc","smc2"),each = nobservations)))
results$time = rep(1:nobservations, 2)
results$logevidence = c(smc_results$logevidence,smc2_results$logevidence)
ggplot(results) + geom_line(aes(time, logevidence/time, color = from), size = 1)

#--------------------------------------------------------------------------------------------
# Check the h-score (RESCALED BY 1/t)
results$hscore = c(smc_results$hscore,smc2_results$hscore)
ggplot(results) + geom_line(aes(time, hscore/time, color = from), size = 1)
