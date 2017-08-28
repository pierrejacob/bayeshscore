#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# Kangaroo Model 1 - Logistic Model
# >> with parameter: # theta = (b, sigma, tau, r)
# >> with prior: independent uniforms
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
model.kangarooLogistic = list()

# dimension of parameter, observations, and states
model.kangarooLogistic$dimtheta = 4
model.kangarooLogistic$dimY = 2
model.kangarooLogistic$dimX = 1

# sampler from the prior distribution on parameters
model.kangarooLogistic$rprior = function(Ntheta){
  b = runif(Ntheta,0,10)
  sigma = runif(Ntheta,0,10)
  tau = runif(Ntheta,0,10)
  r = runif(Ntheta,-10,10)
  return (cbind(b,sigma,tau,r))
}

# density the prior distribution on parameters
model.kangarooLogistic$dprior = function(theta){
  b = theta[1]
  sigma = theta[2]
  tau = theta[3] 
  r = theta[4]
  return (dunif(b,0,10)*dunif(sigma,0,10)*dunif(tau,0,10)*dunif(r,-10,10))
}

# sampler from the initial distribution of the states
model.kangarooLogistic$rinitial = function(theta,N=1){
  return (rlnorm(N,meanlog = 5,sdlog = sqrt(10)))
}

#WARNING: single trajectory !!! (not vectorized)
model.kangarooLogistic$rtransition = function(Xt,dt,theta){
  b = theta[1]
  sigma = theta[2]
  r = theta[4]
  Xold = Xt
  remaining = dt
  while (remaining > 0) {
    A = b*Xold-r-(sigma^2)/2
    B = 6.4*sigma
    if (A>0) {
      delta = min(remaining,((-B+sqrt(B^2+4*A))/(2*A))^2) #max stepsize for Xnew > 0 with high proba
    }
    else {
      delta = min(remaining,((-B+sqrt(B^2-4*A))/(-2*A))^2) #max stepsize for Xnew > 0 with high proba
    }
    Xnew = Xold*(1+(r+(sigma^2)/2-b*Xold)*delta + sigma*sqrt(delta)*rnorm(1,mean = 0,sd = 1))
    Xold = Xnew
    remaining = remaining-delta
  }
  return (Xnew)
}

# density of the observations
model.kangarooLogistic$dobs = function(Yt,Xt,t,theta){
  tau = theta[3]
  n = 1/tau
  p = 1/(1+tau*Xt)
  return (dnbinom(Yt[1],size = n,prob = p)*dnbinom(Yt[2],size = n,prob = p))
}

# OPTIONAL: simulate observations
model.kangarooLogistic$robs = function(Xt,t,theta){
  tau = theta[3]
  n = 1/tau
  p = 1/(1+tau*Xt)
  return (rnbinom(2,size = n,prob = p))
}

# OPTIONAL: lower and upper bounds of observations
model.kangarooLogistic$a = c(0,0)
model.kangarooLogistic$b = c(Inf,Inf)

