#'@rdname model.kangarooExponential
#'@title model.kangarooExponential
#'@description This implements the exponential model in Knape et al. (2012)
#'@export
get_model_kangarooExponential <- function(timesteps = data_kangaroo["time",],range_r = 10){
  model.kangarooExponential = list()
  # dimension of parameter, observations, and states
  model.kangarooExponential$dimtheta = 3
  model.kangarooExponential$dimY = 2
  model.kangarooExponential$dimX = 1
  # sampler from the prior distribution on parameters
  model.kangarooExponential$rprior = function(Ntheta){
    sigma = runif(Ntheta,0,10)
    tau = runif(Ntheta,0,10)
    r = runif(Ntheta,-range_r,range_r)
    return (cbind(sigma,tau,r))
  }
  # density the prior distribution on parameters
  model.kangarooExponential$dprior = function(theta, log = TRUE){
    sigma = theta[1]
    tau = theta[2]
    r = theta[3]
    if (log==TRUE){
      return (dunif(sigma,0,10,log) + dunif(tau,0,10,log) + dunif(r,-range_r,range_r,log))
    }
    else{
      return (dunif(sigma,0,10,log) * dunif(tau,0,10,log) * dunif(r,-range_r,range_r,log))
    }
  }
  # sampler from the initial distribution of the states
  model.kangarooExponential$rinitial = function(theta,N){
    return (matrix(rlnorm(N,meanlog = 5,sdlog = sqrt(10)),nrow = N))
  }
  # sampler from the transition density of the states
  model.kangarooExponential$rtransition = function(Xt,t,theta){
    sigma = theta[1]
    r = theta[3]
    N = nrow(Xt)
    dt = timesteps[t] - timesteps[t-1]
    return (matrix(Xt*exp(r*dt+sigma*rnorm(N,mean = 0,sd = sqrt(dt))),nrow = N))
  }
  # density of the observations
  model.kangarooExponential$dobs = function(Yt,Xt,t,theta,log = TRUE){
    tau = theta[2]
    n = 1/tau
    if (log==TRUE){
      return (dnbinom(Yt[1],size = n,mu = Xt,log = log) + dnbinom(Yt[2],size = n,mu = Xt,log = log))
    }
    else{
      return (dnbinom(Yt[1],size = n,mu = Xt,log = log) * dnbinom(Yt[2],size = n,mu = Xt,log = log))
    }
  }
  # OPTIONAL: simulate observations
  model.kangarooExponential$robs = function(Xt,t,theta){
    tau = theta[2]
    n = 1/tau
    return (rnbinom(2,size = n,mu = Xt))
  }
  # OPTIONAL: lower and upper bounds of observations
  model.kangarooExponential$lower = c(0,0)
  model.kangarooExponential$upper = c(Inf,Inf)
  return(model.kangarooExponential)
}
