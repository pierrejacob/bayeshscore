# compare Rcpp implementtation of Levy-driven SV transition with naive implementation
rm(list = ls())
library(HyvarinenSSM)
library(microbenchmark)
set.seed(19)

nobservations = 10
timesteps = 1:nobservations
theta = c(0, 0, 0.5, 0.0625, 0.01)
model = get_model_SVLevy_singlefactor(timesteps)
modelnaive = model

# Generate some observations
observations = simulateData(model,theta,nobservations)
Y = observations$Y
Nx = 2^10

#---------------------------------------------------------------------------------------
# Test initialization
#---------------------------------------------------------------------------------------
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
microbenchmark(model$rinitial(theta,Nx),modelnaive$rinitial(theta,Nx))

#---------------------------------------------------------------------------------------
# Test transition
#---------------------------------------------------------------------------------------
modelnaive$rtransition = function(Xs,t,theta){
  # parse parameters
  xi = theta[3]
  w2 = theta[4]
  lambda = theta[5]
  Nx = ncol(Xs)
  # random number of jumps, jump times, jump sizes
  k = rpois(Nx,lambda*(xi^2)/w2)
  c1k = lapply(k, function(nbjump) runif(nbjump,timesteps[t-1],timesteps[t]))
  e1k = lapply(k, function(nbjump) rexp(nbjump,xi/w2))
  # propagate
  zt = sapply(1:Nx,function(i){exp(-lambda)*Xs[2,i] + sum(exp(-lambda*(timesteps[t]-c1k[[i]]))*e1k[[i]])})
  vt = sapply(1:Nx, function(i){(1/lambda)*(Xs[2,i] - zt[i] + sum(e1k[[i]]))})
  return (rbind(vt,zt))
}
Xts = model$rinitial(theta,Nx)
microbenchmark(model$rtransition(Xts,2,theta),modelnaive$rtransition(Xts,2,theta))

#---------------------------------------------------------------------------------------
# Test density evaluation
#---------------------------------------------------------------------------------------
modelnaive$dobs = function(Yt,Xts,t,theta,log = TRUE){
  mu = theta[1]
  beta = theta[2]
  Nx = ncol(Xts)
  return (sapply(1:Nx,function(i){dnorm(Yt, mu + beta*Xts[1,i], sqrt(Xts[1,i]), log)}))
}
Xts = model$rinitial(theta,Nx)
microbenchmark(model$dobs(Y[,1],Xts,1,theta),modelnaive$dobs(Y[,1],Xts,1,theta))

#---------------------------------------------------------------------------------------
# Test derivatives evaluation
#---------------------------------------------------------------------------------------
modelnaive$derivativelogdobs = function(Yt,Xts,t,theta){
  mu = theta[1]
  beta = theta[2]
  Nx = ncol(Xts)
  d1 = matrix(sapply(1:Nx,function(i){(mu+beta*Xts[1,i] - Yt)/Xts[1,i]}), ncol = model$dimY)
  d2 = matrix(sapply(1:Nx,function(i){-1/Xts[1,i]}), ncol = model$dimY)
  return (list(jacobian = d1, hessiandiag = d2))
}
Xts = model$rinitial(theta,Nx)
microbenchmark(model$derivativelogdobs(Y[,1],Xts,1,theta),modelnaive$derivativelogdobs(Y[,1],Xts,1,theta))
