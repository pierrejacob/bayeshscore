##################################################################################################
# This checks that the outputs using SMC and SMC2 match the exact results
# in a linear gaussian case with 1 parameter (with discrete prior so that everything can be
# computed exactly).
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(29)

# Define model and data
nobservations = 30
model = get_model_lineargaussiansimpler_discreteprior()
true_sigmav2 = 1
sim = simulateData(model, theta = true_sigmav2, nobservations)
X = sim$X
Y = sim$Y
observations = matrix(Y, nrow = model$dimY)# observations in a matrix of dimensions dimy x nobservations

#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$Nx = 2^7
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$store_X = FALSE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$min_acceptance_rate = 0.0
algorithmic_parameters$nmoves = 1
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R

#--------------------------------------------------------------------------------------------
### Run SMC
smc_results = hscore(observations, model, algorithmic_parameters)
### Run SMC_2
module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree
model_withoutlikelihood = model
model_withoutlikelihood$likelihood = NULL
model_withoutlikelihood$dpredictive = NULL
smc2_results = hscore(observations, model_withoutlikelihood, algorithmic_parameters)

#--------------------------------------------------------------------------------------------
#########################################################################################
#########################################################################################
########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################
thetas_smc = smc_results$thetas_history[[nobservations+1]]
normw_smc = smc_results$normw_history[[nobservations+1]]
#
thetas_smc2 = smc2_results$thetas_history[[nobservations+1]]
normw_smc2 = smc2_results$normw_history[[nobservations+1]]

#--------------------------------------------------------------------------------------------
#Compute exact posterior
phi = model$phi
psi = model$psi
sigmaW2 = model$sigmaW2
initial_mean = 0
initial_var = (sigmaW2)/(1-phi^2)
X_t_t_1 = matrix(NA,nrow = length(model$supportprior),ncol = nobservations)
P_t_t_1 = matrix(NA,nrow = length(model$supportprior),ncol = nobservations)
for (i in 1:length(model$supportprior)){
  kalman = KF_filtering(Y,phi,psi,model$supportprior[i],sigmaW2,initial_mean,initial_var)
  X_t_t_1[i,] = sapply(1:nobservations,function(t)kalman[[t]]$muX_t_t_1)
  P_t_t_1[i,] = sapply(1:nobservations,function(t)kalman[[t]]$PX_t_t_1)
}
posterior_exact = matrix(NA,nrow = length(model$supportprior), ncol = nobservations)
for (i in 1:length(model$supportprior)){
  sigV2 =  model$supportprior[i]
  posterior_exact[i,1] = dnorm(Y[,1],psi*X_t_t_1[i,1],sqrt(sigV2 + (psi^2)*P_t_t_1[i,1]))
}
posterior_exact[,1] = (posterior_exact[,1])/sum(posterior_exact[,1])
for (t in 2:nobservations){
  for (i in 1:length(model$supportprior)){
    sigV2 =  model$supportprior[i]
    posterior_exact[i,t] = posterior_exact[i,t-1]*dnorm(Y[,t],psi*X_t_t_1[i,t],sqrt(sigV2 + (psi^2)*P_t_t_1[i,t]))
  }
  posterior_exact[,t] = (posterior_exact[,t])/sum(posterior_exact[,t])
}
# Checking sample from the posterior distribution (marginal histogram)
Ntheta = algorithmic_parameters$Ntheta
post = data.frame(from = factor(rep(c("smc","smc2"),each = Ntheta)))
post$theta = c(thetas_smc[1,],thetas_smc2[1,])
post$weight = c(normw_smc,normw_smc2)
dx = c(model$supportprior[1],diff(model$supportprior)) #used to renormalize discrete exact posterior (by approximating area by riemann sum)
ggplot(post) +  geom_density(aes(theta, weight = weight, fill = from), alpha = 0.6) +
  geom_line(aes(x,y/dx),data = data.frame(x=model$supportprior,y=posterior_exact[,nobservations]),colour="blue",size=1.5,linetype=1)

#--------------------------------------------------------------------------------------------
#Compute exact log-predictive and evidence
log_py_t_t_1 = rep(NA,nobservations)
log_py_t_t_1_func = list()
log_py_t_t_1_func[[1]] = function(y){
  temp = 0
  for (i in 1:length(model$supportprior)){
    sigV2 =  model$supportprior[i]
    temp = temp + (1/length(model$supportprior))*dnorm(y,phi*X_t_t_1[i,1],sqrt(sigV2 + (psi^2)*P_t_t_1[i,1]))
  }
  return (log(temp))
}
log_py_t_t_1[1] = log_py_t_t_1_func[[1]](Y[,1])
for (t in 2:nobservations){
  log_py_t_t_1_func[[t]] = function(y){
    temp = 0
    for (i in 1:length(model$supportprior)){
      sigV2 =  model$supportprior[i]
      temp = temp + posterior_exact[i,t-1]*dnorm(y,phi*X_t_t_1[i,t],sqrt(sigV2 + (psi^2)*P_t_t_1[i,t]))
    }
    return (log(temp))
  }
  log_py_t_t_1[t] = log_py_t_t_1_func[[t]](Y[,t])
}
logevidence_exact = cumsum(log_py_t_t_1)

# Check the log-evidence (RESCALED BY 1/t)
results = data.frame(from = factor(rep(c("smc","smc2"),each = nobservations)))
results$time = rep(1:nobservations, 2)
results$logevidence = c(smc_results$logevidence,smc2_results$logevidence)
ggplot() +
  geom_line(aes(1:nobservations, logevidence_exact/(1:nobservations)),color="blue",size=2,linetype=2) +
  geom_line(data = results,aes(time, logevidence/time, color = from), size = 1)

#--------------------------------------------------------------------------------------------
#Compute exact (numerically via NumDeriv) hscore
incr_hscore = rep(NA,nobservations)
for (t in 1:nobservations){
  incr_hscore[t] = 2*hessian(log_py_t_t_1_func[[t]],Y[,t]) + (grad(log_py_t_t_1_func[[t]],Y[,t]))^2
}
hscore_exact = cumsum(incr_hscore)

# Check the h-score (RESCALED BY 1/t)
results$hscore = c(smc_results$hscore,smc2_results$hscore)
ggplot() +
  geom_line(aes(1:nobservations, hscore_exact/(1:nobservations)),color="blue",size=2,linetype=2) +
  geom_line(data = results,aes(time, hscore/time, color = from), size = 1)

