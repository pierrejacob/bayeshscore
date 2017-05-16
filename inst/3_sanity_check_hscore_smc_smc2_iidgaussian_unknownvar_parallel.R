rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
set.seed(29)

# Define model and data
nobservations <- 30
model <- get_model_iid_gaussian_unknown_variance()
true_sigmav2 = 1
sim = simulateData(model, theta = true_sigmav2, nobservations)
X = sim$X
Y = sim$Y
observations <- matrix(Y, nrow = model$dimY)
# observations in a matrix of dimensions dimy x nobservations

#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters <- list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$Nx = 2^0
algorithmic_parameters$verbose = FALSE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$store_X = FALSE
algorithmic_parameters$ess_threshold = 0.9 # purposely set high to force some resample-move steps
algorithmic_parameters$min_acceptance_rate = 0.5
algorithmic_parameters$nmoves = 2
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R

#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------
### Run SMC
smc_results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
  hscore(observations, model, algorithmic_parameters)
}
# remove likelihood to force the use of SMC2
model_withoutlikelihood = model
model_withoutlikelihood$likelihood = NULL
model_withoutlikelihood$dpredictive = NULL
### Run SMC2
smc2_results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
  module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
  TreeClass <<- module_tree$Tree
  hscore(observations, model_withoutlikelihood, algorithmic_parameters)
}

#--------------------------------------------------------------------------------------------
#########################################################################################
#########################################################################################
########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################
#########################################################################################
#########################################################################################
#--------------------------------------------------------------------------------------------
#Compute exact posterior
nu = model$df
nu_post = nu + nobservations
s2_post = (nu+sum(Y^2))/nu_post
# Checking sample from the posterior distribution (marginal histogram)
Ntheta = algorithmic_parameters$Ntheta
post = data.frame(from = factor(rep(c("smc","smc2"),each = repl*Ntheta)))
post$theta = c(c(sapply(1:repl,function(i)smc_results[[i]]$thetas_history[[nobservations+1]])),
               c(sapply(1:repl,function(i)smc2_results[[i]]$thetas_history[[nobservations+1]])))
post$weight = c(c(sapply(1:repl,function(i)smc_results[[i]]$normw_history[[nobservations+1]])),
                c(sapply(1:repl,function(i)smc2_results[[i]]$normw_history[[nobservations+1]])))
post$repl = rep(rep(1:repl,each = Ntheta),2)
ggplot(post) +
  geom_density(aes(theta, weight = weight, fill = from, group = interaction(repl,from)), alpha = 0.6) +
  stat_function(fun = function(y)dinvchisq(y,nu_post,s2_post,FALSE),colour="blue",size=1.5,linetype=1)

#--------------------------------------------------------------------------------------------
#compute exact log-evidence
logevidence_exact = rep(NA,nobservations)
for (t in 1:nobservations) {
  nu_t = nu + (t-1)
  st2 = (nu + sum(Y[,1:(t-1)]^2))/nu_t
  logevidence_exact[t] = dtscaled(Y[,t],nu_t,st2,TRUE)
}
logevidence_exact = cumsum(logevidence_exact)
# Check the log-evidence (RESCALED BY 1/t)
results = data.frame(from = factor(rep(c("smc","smc2"),each = repl*nobservations)))
results$repl = rep(rep(1:repl,each = nobservations),2)
results$time = rep(1:nobservations, repl*2)
results$logevidence = c(sapply(1:repl,function(i)smc_results[[i]]$logevidence),
                        sapply(1:repl,function(i)smc2_results[[i]]$logevidence))
ggplot() +
  geom_line(aes(1:nobservations, logevidence_exact/(1:nobservations)),color="blue",size=2,linetype=2) +
  geom_line(data = results,aes(time, logevidence/time, color = from, group = interaction(repl,from)))


#--------------------------------------------------------------------------------------------
#Compute exact h-score
hscore_exact = rep(NA,nobservations)
for (t in 1:nobservations){
  s = sum(Y[,1:t]^2)
  hscore_exact[t] = ((nu+t)/((nu+s)^2))*((nu+t+4)*Y[,t]^2-2*(nu+s))
}
hscore_exact = cumsum(hscore_exact)
# Check the h-score (RESCALED BY 1/t)
results$hscore = c(sapply(1:repl,function(i)smc_results[[i]]$hscore),
                   sapply(1:repl,function(i)smc2_results[[i]]$hscore))
ggplot() +
  geom_line(aes(1:nobservations, hscore_exact/(1:nobservations)),color="blue",size=2,linetype=2) +
  geom_line(data = results,aes(time, hscore/time, color = from, group = interaction(repl,from)))

