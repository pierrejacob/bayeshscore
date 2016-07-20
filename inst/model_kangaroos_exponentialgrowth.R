#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# Kangaroo Model 2 - Exponential Growth
# >> with parameter: # theta = (sigma, tau, r)
# >> with prior: independent uniforms
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
model.kangarooExponential = list()

# dimension of parameter, observations, and states
model.kangarooExponential$dimtheta = 3
model.kangarooExponential$dimY = 2
model.kangarooExponential$dimX = 1


# sampler from the prior distribution on parameters
model.kangarooExponential$rprior = function(Ntheta){
  sigma = runif(Ntheta,0,10)
  tau = runif(Ntheta,0,10)
  r = runif(Ntheta,-10,10)
  return (cbind(sigma,tau,r))
}

# density the prior distribution on parameters
model.kangarooExponential$dprior = function(theta){
  sigma = theta[1]
  tau = theta[2] 
  r = theta[3]
  return (dunif(sigma,0,10)*dunif(tau,0,10)*dunif(r,-10,10))
}

# sampler from the initial distribution of the states
model.kangarooExponential$rinitial = function(theta,N=1){
  return (rlnorm(N,meanlog = 5,sdlog = sqrt(10)))
}

# sampler from the transition density of the states
model.kangarooExponential$rtransition = function(Xt,dt,theta){
  sigma = theta[1]
  r = theta[3]
  N = length(Xt)
  return (Xt*exp(r*dt+sigma*rnorm(N,mean = 0,sd = sqrt(dt))))
}

# density of the observations
model.kangarooExponential$dobs = function(Yt,Xt,t,theta){
  tau = theta[2]
  n = 1/tau
  p = 1/(1+tau*Xt)
  return (dnbinom(Yt[1],size = n,prob = p)*dnbinom(Yt[2],size = n,prob = p))
}

# OPTIONAL: simulate observations
model.kangarooExponential$robs = function(Xt,t,theta){
  tau = theta[2]
  n = 1/tau
  p = 1/(1+tau*Xt)
  return (rnbinom(2,size = n,prob = p))
}

# OPTIONAL: lower and upper bounds of observations
model.kangarooExponential$a = c(0,0)
model.kangarooExponential$b = c(Inf,Inf)


