##################################################################################################
# This checks that the outputs using SMC and SMC2 match the exact results
# in an iid Normal case (with conjugate prior so that everything can be computed analytically).
# Both cases are replicated 'repl' number of times in parallel.
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
set.seed(29)

# Define model and data
nobservations = 20
nu0 = 1
sigma02 = 1
true_sigmav2 = 1
model = get_model_iid_gaussian_unknown_variance(nu0,sigma02)
observations = simulateData(model, theta = true_sigmav2, nobservations)
# observations in a matrix of dimensions dimy x nobservations


#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^7
algorithmic_parameters$verbose = TRUE

algorithmic_parameters$reduce_variance = TRUE
algorithmic_parameters$Nc = 2^10
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R

#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------
### Run SMC
smc_results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
  hscore(observations, get_model_iid_gaussian_unknown_variance(nu0,sigma02), algorithmic_parameters)
}
#--------------------------------------------------------------------------------------------
#########################################################################################
#########################################################################################
########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################
#########################################################################################
#########################################################################################
#--------------------------------------------------------------------------------------------
# #Compute exact posterior
# nu_post = nu0 + nobservations
# s2_post = (nu0+sum(observations^2))/nu_post
# # Checking sample from the posterior distribution (marginal histogram)
# Ntheta = algorithmic_parameters$Ntheta
# post = data.frame(from = factor(rep(c("smc"),each = repl*Ntheta)))
# post$theta = c(c(sapply(1:repl,function(i)smc_results[[i]]$thetas)))
# post$weight = c(c(sapply(1:repl,function(i)smc_results[[i]]$normw)))
# post$repl = rep(1:repl,each = Ntheta)
# ggplot(post) +
#   geom_density(aes(theta, weight = weight, fill = from, group = interaction(repl,from)), alpha = 0.6) +
#   stat_function(fun = function(y)dinvchisq(y,nu_post,s2_post,FALSE),colour="blue",size=1.5,linetype=1)
#
# #--------------------------------------------------------------------------------------------
# #compute exact log-evidence
# logevidence_exact = rep(NA,nobservations)
# for (t in 1:nobservations) {
#   nu_t = nu0 + (t-1)
#   st2 = (nu0 + sum(observations[,1:(t-1)]^2))/nu_t
#   logevidence_exact[t] = dtscaled(observations[,t],nu_t,st2,TRUE)
# }
# logevidence_exact = cumsum(logevidence_exact)
# # Check the log-evidence (RESCALED BY 1/t)
results = data.frame(from = factor(rep(c("smc"),each = repl*nobservations)))
results$repl = rep(1:repl,each = nobservations)
results$time = 1:nobservations
results$logevidence = c(sapply(1:repl,function(i)smc_results[[i]]$logevidence))
# ggplot() +
#   geom_line(aes(1:nobservations, logevidence_exact/(1:nobservations)),color="blue",size=2,linetype=2) +
#   geom_line(data = results,aes(time, logevidence/time, color = from, group = interaction(repl,from)))
#

#--------------------------------------------------------------------------------------------
#Compute exact h-score
hscore_exact = rep(NA,nobservations)
for (t in 1:nobservations){
  s = sum(observations[,1:t]^2)
  hscore_exact[t] = ((nu0+t)/((nu0+s)^2))*((nu0+t+4)*observations[,t]^2-2*(nu0+s))
}
hscore_exact = cumsum(hscore_exact)
# Check the h-score (RESCALED BY 1/t)
results$hscore = c(sapply(1:repl,function(i)smc_results[[i]]$hscore))
ggplot() +
  geom_line(aes(1:nobservations, hscore_exact/(1:nobservations)),color="blue",size=2,linetype=2) +
  geom_line(data = results,aes(time, hscore/time, color = from, group = repl))

