rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(29)

# create data
nobservations <- 50
model <- get_model_simplerlineargaussian()
theta_star <- c(0.8,1)
sim = simulateData(model, theta = theta_star, nobservations)
X = sim$X
Y = sim$Y
observations <- matrix(Y, nrow = model$dimY) # observations in a matrix of dimensions dimy x nobservations

# set algorithmic parameters
algorithmic_parameters <- list()
algorithmic_parameters$Ntheta = 1024
algorithmic_parameters$Nx = 128
algorithmic_parameters$observation_type = 'continuous'
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$min_acceptance_rate = 0.45
algorithmic_parameters$nmoves = 2

###
model$likelihood = NULL #forces the use of SMC2
model_numderiv = model
model_numderiv$derivativelogdobs = NULL #forces the use of numerical derivative
###

module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree
smc2_results <- hscore(observations, model, algorithmic_parameters)
smc2_results_numderiv <- hscore(observations, model_numderiv, algorithmic_parameters)

time_t = 50
thetas_smc2 <- smc2_results$thetas_history[[time_t+1]]
normw_smc2 <- smc2_results$normw_history[[time_t+1]]
thetas_smc2numderiv <- smc2_results_numderiv$thetas_history[[time_t+1]]
normw_smc2numderiv <- smc2_results_numderiv$normw_history[[time_t+1]]
# Checking sample from the posterior distribution (marginal histogram)
Ntheta = algorithmic_parameters$Ntheta
post = data.frame(from = factor(rep(c("smc2","smc2numderiv"),each = Ntheta)))
post$theta1 = c(thetas_smc2[1,],thetas_smc2numderiv[1,])
post$theta2 = c(thetas_smc2[2,],thetas_smc2numderiv[2,])
post$weight = c(normw_smc2,normw_smc2numderiv)
plot_theta1 = ggplot(post) +  geom_density(aes(theta1, weight = weight, fill = from), alpha = 0.6) + theme(legend.position="none")
plot_theta2 = ggplot(post) +  geom_density(aes(theta2, weight = weight, fill = from), alpha = 0.6)
grid.arrange(plot_theta1, plot_theta2, ncol = 2, widths=c(1,1.35))
# Check the log-evidence (RESCALED BY 1/t)
results = data.frame(from = factor(rep(c("smc2","smc2numderiv"),each = time_t)))
results$time = rep(1:time_t, 2)
results$logevidence = c(smc2_results$logevidence,smc2_results_numderiv$logevidence)
ggplot(results) + geom_line(aes(time, logevidence/time, color = from), size = 1)
# Check the h-score (RESCALED BY 1/t)
ggplot() +
  geom_line(aes(1:time_t,smc2_results$hscore[1:time_t]/1:time_t), color = "purple", size = 1) +
  geom_line(aes(1:time_t,smc2_results_numderiv$hscore[1:time_t]/1:time_t), color = "blue", size = 1)
