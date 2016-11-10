#'@rdname get_model_kangarooRandomwalk
#'@title get_model_kangarooRandomwalk
#'@description This implements the random walk model in Knape et al. (2012)
#'@export
get_model_kangarooRandomwalk <- function(timesteps = data_kangaroo["time",]){
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
  model.kangarooRandomwalk$dprior = function(theta, log = TRUE){
    sigma = theta[1]
    tau = theta[2]
    if (log==TRUE){
      return (dunif(sigma,0,10,log) + dunif(tau,0,10,log))
    }
    else{
      return (dunif(sigma,0,10,log) * dunif(tau,0,10,log))
    }
  }
  # sampler from the initial distribution of the states
  model.kangarooRandomwalk$rinitial = function(theta,N){
    return (matrix(rlnorm(N,meanlog = 5,sdlog = sqrt(10)), nrow = N))
  }
  # sampler from the transition density of the states
  model.kangarooRandomwalk$rtransition = function(Xt,t,theta){
    sigma = theta[1]
    N = nrow(Xt)
    dt = timesteps[t] - timesteps[t-1]
    return (matrix(Xt*exp(sigma*rnorm(N,mean = 0,sd = sqrt(dt))), nrow = N))
  }
  # density of the observations
  model.kangarooRandomwalk$dobs = function(Yt,Xt,t,theta,log = TRUE){
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
  model.kangarooRandomwalk$robs = function(Xt,t,theta){
    tau = theta[2]
    n = 1/tau
    return (rnbinom(2,size = n,mu = Xt))
  }
  # OPTIONAL: lower and upper bounds of observations
  model.kangarooRandomwalk$lower = c(0,0)
  model.kangarooRandomwalk$upper = c(Inf,Inf)
  return(model.kangarooRandomwalk)
}
