#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# Kangaroo Model 3 - log-Random Walk
# >> with parameter: # theta = (sigma, tau)
# >> with prior: independent uniforms
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
model.kangarooRandomwalk = list()

# dimension of parameter, observations, and states
model.kangarooRandomwalk$dimtheta = 2
model.kangarooRandomwalk$dimY = 2
model.kangarooRandomwalk$dimX = 1

# sampler from the prior distribution on parameters
model.kangarooRandomwalk$rprior = function(Ntheta){
  sigma = runif(Ntheta,0,10)
  tau = runif(Ntheta,0,10)
  return (cbind(sigma,tau))
}

# density the prior distribution on parameters
model.kangarooRandomwalk$dprior = function(theta){
  sigma = theta[1]
  tau = theta[2]
  return (dunif(sigma,0,10)*dunif(tau,0,10))
}

# sampler from the initial distribution of the states
model.kangarooRandomwalk$rinitial = function(theta,N=1){
  return (rlnorm(N,meanlog = 5,sdlog = sqrt(10)))
}

# sampler from the transition density of the states
model.kangarooRandomwalk$rtransition = function(Xt,dt,theta){
  sigma = theta[1]
  N = length(Xt)
  return (Xt*exp(sigma*rnorm(N,mean = 0,sd = sqrt(dt))))
}

# density of the observations
model.kangarooRandomwalk$dobs = function(Yt,Xt,t,theta){
  tau = theta[2]
  n = 1/tau
  p = 1/(1+tau*Xt)
  return (dnbinom(Yt[1],size = n,prob = p)*dnbinom(Yt[2],size = n,prob = p))
}

# OPTIONAL: simulate observations
model.kangarooRandomwalk$robs = function(Xt,t,theta){
  tau = theta[2]
  n = 1/tau
  p = 1/(1+tau*Xt)
  return (rnbinom(2,size = n,prob = p))
}

# OPTIONAL: lower and upper bounds of observations
model.kangarooRandomwalk$a = c(0,0)
model.kangarooRandomwalk$b = c(Inf,Inf)
