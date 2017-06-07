library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(wesanderson)
set.seed(19)
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$Nx = 2^3
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$save = FALSE
algorithmic_parameters$store_X_history = TRUE
algorithmic_parameters$store_thetas_history = TRUE
algorithmic_parameters$store_last_X = FALSE
algorithmic_parameters$store_last_thetas = TRUE


nobservations = 3
timesteps = 1:nobservations
theta = c(0, 0, 0.5, 0.0625, 0.01)
# observations = simulateData(get_model_SVLevy_singlefactor(timesteps),theta,nobservations)$Y
observations = matrix(c(-0.1549653,0.08272682,-0.1054415), nrow = 1)


nobservations = ncol(observations)
timesteps = 1:nobservations
# define models
model = get_model_SVLevy_singlefactor(timesteps)
results = hscore(observations, model, algorithmic_parameters)


mu = -0.09948713
beta = 3.38691535
xi = 0.00602995
w2 = 3.54957726
lambda = 0.05853343

hist(model$rprior(algorithmic_parameters$Ntheta)[1,],breaks=100)
hist(model$rprior(algorithmic_parameters$Ntheta)[2,],breaks=100)
hist(model$rprior(algorithmic_parameters$Ntheta)[3,],breaks=100)
hist(model$rprior(algorithmic_parameters$Ntheta)[4,],breaks=100)
hist(model$rprior(algorithmic_parameters$Ntheta)[5,],breaks=100)

theta = c(mu, beta, xi, w2, lambda)
model$rinitial(theta,algorithmic_parameters$Nx)

hist(model$rinitial(theta,2^10)[1,],breaks=100)
hist(model$rinitial(theta,2^10)[2,],breaks=100)

modelnaive = model
modelnaive$rinitial = function(theta,Nx){
  # parse parameters
  xi = theta[3]
  w2 = theta[4]
  lambda = theta[5]
  # initial z0
  z0 = rgamma(Nx,shape = (xi^2)/w2,rate = xi/w2)
  # random number of jumps, jump times, jump sizes
  k = rpois(Nx,lambda*(xi^2)/w2)
  c1k = lapply(k, function(nbjump) runif(nbjump,0,timesteps[1]))
  e1k = lapply(k, function(nbjump) rexp(nbjump,xi/w2))
  # propagate
  z1 = sapply(1:Nx,function(i){exp(-lambda)*z0[i] + sum(exp(-lambda*(timesteps[1]-c1k[[i]]))*e1k[[i]])})
  v1 = sapply(1:Nx, function(i){(1/lambda)*(z0[i] - z1[i] + sum(e1k[[i]]))})
  return (rbind(v1,z1))
}
hist(modelnaive$rinitial(theta,2^10)[1,],breaks=100)
hist(modelnaive$rinitial(theta,2^10)[2,],breaks=100)

Nx = 1000
sum(rpois(Nx,lambda*(xi^2)/w2))
sum(rgamma(Nx,shape = (xi^2)/w2, rate = xi/w2))


