rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(19)

# Define model and data
nobservations <- 5
model <- get_model_poissonhmm()
true_theta = 0.3
sim = simulateData(model, theta = true_theta, nobservations)
X = sim$X
Y = sim$Y
observations <- matrix(Y, nrow = model$dimY)
# observations in a matrix of dimensions dimy x nobservations

#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters <- list()
algorithmic_parameters$Ntheta = 2^6
algorithmic_parameters$Nx = 2^7
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$store_X = FALSE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$min_acceptance_rate = 0.3
algorithmic_parameters$nmoves = 1
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R

#--------------------------------------------------------------------------------------------
### Run SMC
#NB: the exact likelihood has been implemented in the most naive (hence inefficient) way
# this code is just meant as a proof of concept
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
thetas_smc <- smc_results$thetas_history[[nobservations+1]]
normw_smc <- smc_results$normw_history[[nobservations+1]]
#
thetas_smc2 <- smc2_results$thetas_history[[nobservations+1]]
normw_smc2 <- smc2_results$normw_history[[nobservations+1]]

#Compute exact likelihood
l = model$lambda
lx = model$lambdaX
p_Y_1_t_theta = function(Y_1_t,theta){
  t = ncol(Y_1_t)
  all_paths = expand.grid(rep(list(0:1), t))
  if (t==1){
    temp = 0
    for (i in 1:(2^t)){
      X = all_paths[i,1:t]
      rate = l+X*lx
      temp = temp + prod(exp(-(rate))*(rate^Y_1_t)/factorial(Y_1_t))*(1/2) + 0*theta
      #the "0*theta" term (artificially) allows integrate to recognize this as a function of theta
    }
    return (temp)
  }
  else {
    temp = 0
    for (i in 1:(2^t)){
      X = all_paths[i,1:t]
      rate = l+X*lx
      temp = temp + prod(exp(-(rate))*(rate^Y_1_t)/factorial(Y_1_t))*(1/2)*(theta^sum((1-X[2:t])*(1-X[1:(t-1)])+X[2:t]*X[1:(t-1)]))*((1-theta)^sum(X[2:t]*(1-X[1:(t-1)])+(1-X[2:t])*X[1:(t-1)]))
    }
    return (temp)
  }
}
#Compute exact evidence
p_Y_1_t = function(Y_1_t) {
  likelihood = function(theta) p_Y_1_t_theta(Y_1_t,theta)
  return (integrate(likelihood,0,1)$value)
}
#Compute exact posterior
normalizing = integrate(function(theta) p_Y_1_t_theta(Y,theta),0,1)$value
post_exact = function(theta){
  return ((p_Y_1_t_theta(Y,theta)/normalizing)*(theta<=1)*(theta>=0))
}
# Checking sample from the posterior distribution (marginal histogram)
Ntheta = algorithmic_parameters$Ntheta
post = data.frame(from = factor(rep(c("smc","smc2"),each = Ntheta)))
post$theta = c(thetas_smc[1,],thetas_smc2[1,])
post$weight = c(normw_smc,normw_smc2)
ggplot(post) +  geom_density(aes(theta, weight = weight, fill = from), alpha = 0.6) +
  stat_function(fun = function(y)post_exact(y),colour="blue",size=1.5,linetype=1)


#--------------------------------------------------------------------------------------------
#Compute exact log-evidence
logevidence_exact = rep(NA,nobservations)
for (t in 1:nobservations) {
  logevidence_exact[t] = log(p_Y_1_t(Y[,1:t,drop=FALSE]))
}
# Check the log-evidence (RESCALED BY 1/t)
results = data.frame(from = factor(rep(c("smc","smc2"),each = nobservations)))
results$time = rep(1:nobservations, 2)
results$logevidence = c(smc_results$logevidence,smc2_results$logevidence)
ggplot() +
  geom_line(aes(1:nobservations, logevidence_exact/(1:nobservations)),color="blue",size=2,linetype=2) +
  geom_line(data = results,aes(time, logevidence/time, color = from), size = 1)

#--------------------------------------------------------------------------------------------
#Compute exact predictive
py_t_t_1_func = list()
for (t in 1:nobservations){
  py_t_t_1_func[[t]] = function(Yt){
    if (t==1){
      return (p_Y_1_t(matrix(Yt,ncol=1)))
    }
    else{
      Y_1_t = Y[,1:t,drop=FALSE]
      Y_1_t[,t] = Yt
      return (p_Y_1_t(Y_1_t)/p_Y_1_t(Y_1_t[,1:(t-1),drop=FALSE]))
    }
  }
}
#Compute exact hscore
hincrement = function(k,a,b,d,y,t) {
  ek = rep(0,d)
  ek[k] = 1
  q = py_t_t_1_func[[t]]
  if (y[k]==b[k]) {
    qy = q(y)
    qy_minusek = q(y-ek)
    return (-2*(qy-qy_minusek)/qy_minusek)
  }
  else {
    if (y[k]==a[k]) {
      qy = q(y)
      qy_plusek = q(y+ek)
      return (2*(qy_plusek-qy)/qy + ((qy_plusek-qy)/qy)^2)
    }
    else {
      qy = q(y)
      qy_minusek = q(y-ek)
      qy_plusek = q(y+ek)
      return (2*((qy_plusek-qy)/qy-(qy-qy_minusek)/qy_minusek) + ((qy_plusek-qy)/qy)^2)
    }
  }
}
#Compute exact hscore
incr_hscore = rep(NA,nobservations)
for (t in 1:nobservations){
  for (k in 1:model$dimY){
    incr_hscore[t] = hincrement(k,model$lower,model$upper,model$dimY,Y[,t],t)
  }
}
hscore_exact = cumsum(incr_hscore)
# Check the h-score (RESCALED BY 1/t)
results$hscore = c(smc_results$hscore,smc2_results$hscore)
ggplot() +
  geom_line(aes(1:nobservations, hscore_exact/(1:nobservations)),color="blue",size=2,linetype=2) +
  geom_line(data = results,aes(time, hscore/time, color = from), size = 1)

