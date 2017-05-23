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
nobservations = 30
{
  # The model uses iid observations but we artificially write it as a (degenerate) state-space model
  # to check the validity of the SMC2
  model = list()
  model$observation_type = 'continuous'
  model$df = 20
  model$sigmaW2 = 0
  # dimension of parameter
  model$dimtheta = 1
  model$dimY = 1
  model$dimX = 1
  # sampler from the prior distribution on parameters
  model$rprior = function(Ntheta){
    return (rbind(rinvchisq(Ntheta,model$df,1)))
  }
  # prior distribution density on parameters
  model$dprior = function(theta, log = TRUE){
    sigmaV2 = theta[1]
    return (dinvchisq(sigmaV2,model$df,1,log))
  }
  # sampler from the initial distribution of the states
  model$rinitial = function(theta,N){
    return (matrix(rnorm(N, mean = 0, sd = sqrt(model$sigmaW2)), ncol = N))
  }
  # sampler from the transition density of the states
  model$rtransition = function(Xt,t,theta){
    N = ncol(Xt)
    return (matrix(rnorm(N, mean = 0, sd = sqrt(model$sigmaW2)), ncol = N))
  }
  # density of the observations
  model$dobs = function(Yt,Xt,t,theta,log = TRUE){
    sigmaV2 = theta[1]
    return (dnorm(Yt,mean = Xt,sd = sqrt(sigmaV2), log))
  }
  # first and second partial derivatives of the observation log-density
  # The function is vectorized with respect to the states Xt (dimX by Nx), so that it outputs:
  # >> the jacobian (Nx by dimY matrix: each row is the transpose of the corresponding gradients row-wise)
  # >> the Hessian diagonals (Nx by dimY matrix: each row is the diagonal coeffs of the corresponding Hessian)
  model$derivativelogdobs = function(Yt,Xt,t,theta,dimY){
    sigmaV2 = theta[1]
    N = ncol(Xt)
    d1 = t((Xt-repeat_column(N,Yt))/sigmaV2)
    d2 = matrix(-1/sigmaV2,nrow = N, ncol = dimY)
    return (list(jacobian = d1, hessiandiag = d2))
  }
  # OPTIONAL: one-step predicitve density of the observation at time t given all the past from 1 to (t-1)
  # This relies on some Kalman filter (passed as a byproduct)
  model$dpredictive = function(observations,t,theta,log = TRUE){
    sigmaV2 = theta[1]
    return (dnorm(observations[,t],mean = 0,sd = sqrt(sigmaV2), log))
  }
  # OPTIONAL: simulate observations
  model$robs = function(Xt,t,theta){
    sigmaV2 = theta[1]
    N = ncol(Xt)
    return (matrix(Xt + rnorm(N, mean = 0, sd = sqrt(sigmaV2)),ncol = N))
  }
}
true_sigmav2 = 1
sim = simulateData(model, theta = true_sigmav2, nobservations)
X = sim$X
Y = sim$Y
observations = matrix(Y, nrow = model$dimY)
# observations in a matrix of dimensions dimy x nobservations

#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$Nx = 2^0
algorithmic_parameters$verbose = TRUE
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
  {
    # The model uses iid observations but we artificially write it as a (degenerate) state-space model
    # to check the validity of the SMC2
    model = list()
    model$observation_type = 'continuous'
    model$df = 20
    model$sigmaW2 = 0
    # dimension of parameter
    model$dimtheta = 1
    model$dimY = 1
    model$dimX = 1
    # sampler from the prior distribution on parameters
    model$rprior = function(Ntheta){
      return (rbind(rinvchisq(Ntheta,model$df,1)))
    }
    # prior distribution density on parameters
    model$dprior = function(theta, log = TRUE){
      sigmaV2 = theta[1]
      return (dinvchisq(sigmaV2,model$df,1,log))
    }
    # sampler from the initial distribution of the states
    model$rinitial = function(theta,N){
      return (matrix(rnorm(N, mean = 0, sd = sqrt(model$sigmaW2)), ncol = N))
    }
    # sampler from the transition density of the states
    model$rtransition = function(Xt,t,theta){
      N = ncol(Xt)
      return (matrix(rnorm(N, mean = 0, sd = sqrt(model$sigmaW2)), ncol = N))
    }
    # density of the observations
    model$dobs = function(Yt,Xt,t,theta,log = TRUE){
      sigmaV2 = theta[1]
      return (dnorm(Yt,mean = Xt,sd = sqrt(sigmaV2), log))
    }
    # first and second partial derivatives of the observation log-density
    # The function is vectorized with respect to the states Xt (dimX by Nx), so that it outputs:
    # >> the jacobian (Nx by dimY matrix: each row is the transpose of the corresponding gradients row-wise)
    # >> the Hessian diagonals (Nx by dimY matrix: each row is the diagonal coeffs of the corresponding Hessian)
    model$derivativelogdobs = function(Yt,Xt,t,theta,dimY){
      sigmaV2 = theta[1]
      N = ncol(Xt)
      d1 = t((Xt-repeat_column(N,Yt))/sigmaV2)
      d2 = matrix(-1/sigmaV2,nrow = N, ncol = dimY)
      return (list(jacobian = d1, hessiandiag = d2))
    }
    # OPTIONAL: one-step predicitve density of the observation at time t given all the past from 1 to (t-1)
    # This relies on some Kalman filter (passed as a byproduct)
    model$dpredictive = function(observations,t,theta,log = TRUE){
      sigmaV2 = theta[1]
      return (dnorm(observations[,t],mean = 0,sd = sqrt(sigmaV2), log))
    }
    # OPTIONAL: simulate observations
    model$robs = function(Xt,t,theta){
      sigmaV2 = theta[1]
      N = ncol(Xt)
      return (matrix(Xt + rnorm(N, mean = 0, sd = sqrt(sigmaV2)),ncol = N))
    }
  }
  hscore(observations, model, algorithmic_parameters)
}
# remove likelihood to force the use of SMC2
### Run SMC2
smc2_results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
  {
    # The model uses iid observations but we artificially write it as a (degenerate) state-space model
    # to check the validity of the SMC2
    model = list()
    model$observation_type = 'continuous'
    model$df = 20
    model$sigmaW2 = 0
    # dimension of parameter
    model$dimtheta = 1
    model$dimY = 1
    model$dimX = 1
    # sampler from the prior distribution on parameters
    model$rprior = function(Ntheta){
      return (rbind(rinvchisq(Ntheta,model$df,1)))
    }
    # prior distribution density on parameters
    model$dprior = function(theta, log = TRUE){
      sigmaV2 = theta[1]
      return (dinvchisq(sigmaV2,model$df,1,log))
    }
    # sampler from the initial distribution of the states
    model$rinitial = function(theta,N){
      return (matrix(rnorm(N, mean = 0, sd = sqrt(model$sigmaW2)), ncol = N))
    }
    # sampler from the transition density of the states
    model$rtransition = function(Xt,t,theta){
      N = ncol(Xt)
      return (matrix(rnorm(N, mean = 0, sd = sqrt(model$sigmaW2)), ncol = N))
    }
    # density of the observations
    model$dobs = function(Yt,Xt,t,theta,log = TRUE){
      sigmaV2 = theta[1]
      return (dnorm(Yt,mean = Xt,sd = sqrt(sigmaV2), log))
    }
    # first and second partial derivatives of the observation log-density
    # The function is vectorized with respect to the states Xt (dimX by Nx), so that it outputs:
    # >> the jacobian (Nx by dimY matrix: each row is the transpose of the corresponding gradients row-wise)
    # >> the Hessian diagonals (Nx by dimY matrix: each row is the diagonal coeffs of the corresponding Hessian)
    model$derivativelogdobs = function(Yt,Xt,t,theta,dimY){
      sigmaV2 = theta[1]
      N = ncol(Xt)
      d1 = t((Xt-repeat_column(N,Yt))/sigmaV2)
      d2 = matrix(-1/sigmaV2,nrow = N, ncol = dimY)
      return (list(jacobian = d1, hessiandiag = d2))
    }
    # OPTIONAL: simulate observations
    model$robs = function(Xt,t,theta){
      sigmaV2 = theta[1]
      N = ncol(Xt)
      return (matrix(Xt + rnorm(N, mean = 0, sd = sqrt(sigmaV2)),ncol = N))
    }
  }
  module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
  TreeClass <<- module_tree$Tree
  hscore(observations, model, algorithmic_parameters)
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

