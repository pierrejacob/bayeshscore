##################################################################################################
# This checks that the outputs using SMC and SMC2 match the exact results
# in an iid Normal case (with conjugate prior so that everything can be computed analytically)
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(29)

# Define model and data
nu0 = 1
s02 = 1
nobservations = 5
observations = matrix(rnorm(nobservations), nrow = 1)# observations in a matrix of dimensions dimy x nobservations
model = get_model_iid_gaussian_unknown_variance(nu0,s02)
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^7
algorithmic_parameters$verbose = TRUE

algorithmic_parameters$use_kde = TRUE
algorithmic_parameters$parameters_kde = list(Ny = 2^10)
#--------------------------------------------------------------------------------------------
### Run SMC
smc_results = hscore(observations, model, algorithmic_parameters)
#--------------------------------------------------------------------------------------------
#########################################################################################
#########################################################################################
########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################
thetas_smc = smc_results$thetas
normw_smc = smc_results$normw
#--------------------------------------------------------------------------------------------
#Compute exact posterior
nu_post = nu0 + nobservations
s2_post = (nu0+sum(observations^2))/nu_post
# Checking sample from the posterior distribution (marginal histogram)
Ntheta = algorithmic_parameters$Ntheta
post = data.frame(from = factor(rep("smc",each = Ntheta)))
post$theta = c(thetas_smc[1,])
post$weight = c(normw_smc)
ggplot(post) +  geom_density(aes(theta, weight = weight, fill = from), alpha = 0.6) +
  stat_function(fun = function(y)dinvchisq(y,nu_post,s2_post,FALSE),colour="blue",size=1.5,linetype=1)

#--------------------------------------------------------------------------------------------
#compute exact log-evidence
logevidence_exact = rep(NA,nobservations)
for (t in 1:nobservations) {
  nu_t = nu0 + (t-1)
  st2 = (nu0 + sum(observations[,1:(t-1)]^2))/nu_t
  logevidence_exact[t] = dtscaled(observations[,t],nu_t,st2,TRUE)
}
logevidence_exact = cumsum(logevidence_exact)
# Check the log-evidence (RESCALED BY 1/t)
results = data.frame(from = factor(rep(c("smc"),each = nobservations)))
results$time = 1:nobservations
results$logevidence = c(smc_results$logevidence)
ggplot() +
  geom_line(aes(1:nobservations, logevidence_exact/(1:nobservations)),color="blue",size=2,linetype=2) +
  geom_line(data = results,aes(time, logevidence/time, color = from), size = 1)


#--------------------------------------------------------------------------------------------
#Compute exact h-score
hscore_exact = rep(NA,nobservations)
for (t in 1:nobservations){
  s = sum(observations[,1:t]^2)
  hscore_exact[t] = ((nu0+t)/((nu0+s)^2))*((nu0+t+4)*observations[,t]^2-2*(nu0+s))
}
hscore_exact = cumsum(hscore_exact)
# Check the h-score (RESCALED BY 1/t)
results$hscore = c(smc_results$hscore)
ggplot() +
  geom_line(aes(1:nobservations, hscore_exact/(1:nobservations)),color="blue",size=2,linetype=2) +
  geom_line(data = results,aes(time, hscore/time, color = from), size = 1) +
  geom_line(aes(1:nobservations, smc_results$hscoreKDE/1:nobservations), size = 1.5, color = "black", linetype=4)

