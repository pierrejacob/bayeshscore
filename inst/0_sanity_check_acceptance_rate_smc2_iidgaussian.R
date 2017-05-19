##################################################################################################
# This checks the acceptance rates of the SMC on an i.i.d. gaussian case
# For this case: the SMC uses the analytical predictive distribution (no byproducts needed)
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
set.seed(29)
module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree
#--------------------------------------------------------------------------------------------
# Define model
muprior = 0
sigma2prior = 100
model = get_model_iid_gaussian_unknown_mean(muprior,sigma2prior)
# Generate some simulated data
nobservations = 10
Y = rnorm(nobservations,0,1)
observations = matrix(Y, nrow = 1)# observations in a matrix of dimensions dimy x nobservations
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$Nx = 2^1
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$nmoves = 5 # purposely set high to check the coherence of acceptance rates
#--------------------------------------------------------------------------------------------
# Remove likelihood and predictive in order to force the use of SMC2
model_nolikelihood = model
model_nolikelihood$likelihood = NULL
model_nolikelihood$dpredictive = NULL
# Run SMC2
smc2_result = hscore(observations, model_nolikelihood, algorithmic_parameters)
#--------------------------------------------------------------------------------------------
sigma2post = rep(NA,nobservations)
mu_post = rep(NA,nobservations)
for (t in 1:nobservations){
  # Compute exact posterior
  sigma2post[t] = 1/(1/sigma2prior + t)
  mu_post[t] = ((1/sigma2prior)*muprior + sum(observations[1:t]))*sigma2post[t]
  # # Checking sample from the posterior distribution
  # plot(ggplot(data.frame(theta = smc2_result$thetas_history[[t+1]][1,],
  #                   weight = smc2_result$normw_history[[t+1]])) +
  #        geom_density(aes(theta, weight=weight),fill='blue',alpha=0.3) +
  #        stat_function(fun = function(y)dnorm(y,mu_post[t],sqrt(sigma2post[t])),colour="blue",size=1.5))
}
# Checking last posterior distribution
plot(ggplot(data.frame(theta = smc2_result$thetas_history[[nobservations+1]][1,],
                       weight = smc2_result$normw_history[[nobservations+1]])) +
       geom_density(aes(theta, weight=weight),fill='blue',alpha=0.3) +
       stat_function(fun = function(y)dnorm(y,mu_post[nobservations],sqrt(sigma2post[nobservations])),colour="blue",size=1.5))
#--------------------------------------------------------------------------------------------
# Compute exact logevidence (available analytically as the normalizing constant of the posterior)
exact_ll_at_0 = cumsum(sapply(1:nobservations,function(t)model$dpredictive(observations,t,0)))
exact_logevidence = sapply(1:nobservations,function(t)model$dprior(0) + exact_ll_at_0[t] - dnorm(0,mu_post[t],sqrt(sigma2post[t]),TRUE))
# Check the log-evidence (RESCALED BY 1/t)
ggplot() +
  geom_line(aes(1:nobservations, smc2_result$logevidence),size = 1) +
  geom_line(aes(1:nobservations, exact_logevidence), col ="blue", linetype = "dashed", size = 2)
#--------------------------------------------------------------------------------------------

