#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# Linear Gaussian State-space model (univariate) 
# >> with parameter: # theta = (phi, psi, sigmaV2, sigmaW2)
# >> with prior: independent uniforms
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
model.lineargaussian = list()

# dimension of parameter
model.lineargaussian$dimtheta = 4
model.lineargaussian$dimY = 1
model.lineargaussian$dimX = 1

# sampler from the prior distribution on parameters
model.lineargaussian$rprior = function(Ntheta){
  phi = runif(Ntheta,0.1,0.9)
  psi = runif(Ntheta,0.5,1.5)
  sigmaV2 = runif(Ntheta,0.1,2)
  sigmaW2 = runif(Ntheta,1,10)
  return (cbind(phi,psi,sigmaV2,sigmaW2))
}

# density the prior distribution on parameters
model.lineargaussian$dprior = function(theta){
  phi = theta[1]
  psi = theta[2]
  sigmaV2 = theta[3] 
  sigmaW2 = theta[4]
  return (dunif(phi,0.1,0.9)*dunif(psi,0.5,1.5)*dunif(sigmaV2,0.1,2)*dunif(sigmaW2,1,10))
}

# sampler from the initial distribution of the states
model.lineargaussian$rinitial = function(theta,N=1){
  phi = theta[1]
  sigmaW2 = theta[4]
  return (rnorm(N, mean = 0, sd = sqrt((sigmaW2)/(1-phi^2))))
}

# sampler from the transition density of the states
model.lineargaussian$rtransition = function(Xt,t,theta){
  phi = theta[1]
  sigmaW2 = theta[4]
  N = length(Xt)
  return (phi*Xt + rnorm(N, mean = 0, sd = sqrt(sigmaW2)))
}

# density of the observations
model.lineargaussian$dobs = function(Yt,Xt,t,theta){
  psi = theta[2]
  sigmaV2 = theta[3]
  return (dnorm(Yt,mean = psi*Xt,sd = sqrt(sigmaV2)))
}

# OPTIONAL: simulate observations
model.lineargaussian$robs = function(Xt,t,theta){
  psi = theta[2]
  sigmaV2 = theta[3]
  N = length(Xt)
  return (psi*Xt + rnorm(N, mean = 0, sd = sqrt(sigmaV2)))
}

# OPTIONAL: fixed parameter for simulation
model.lineargaussian$theta = c(0.8,1,4,1)

# OPTIONAL: initial mean and variance for state (for linear gaussian model only)
model.lineargaussian$initialmean = 0
model.lineargaussian$initialvar = (model.lineargaussian$theta[4])/(1-model.lineargaussian$theta[1]^2)
