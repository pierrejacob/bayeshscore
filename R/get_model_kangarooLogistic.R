#'@rdname get_model_kangarooRandomwalk
#'@title get_model_kangarooRandomwalk
#'@description This implements the logistic model in Knape et al. (2012)
#'@export
get_model_kangarooLogistic <- function(){
  model.kangarooLogistic = list()
  # dimension of parameter and variables
  model.kangarooLogistic$dimtheta = 4
  model.kangarooLogistic$dimY = 2
  model.kangarooLogistic$dimX = 1
  # sampler from the prior distribution on parameters
  model.kangarooLogistic$rprior = function(Ntheta){
    sigma = runif(Ntheta,0,10)
    tau = runif(Ntheta,0,10)
    r = runif(Ntheta,-10,10)
    b = runif(Ntheta,0,10)
    return (cbind(sigma,tau,r,b))
  }
  # density the prior distribution on parameters
  model.kangarooLogistic$dprior = function(theta, log = TRUE){
    sigma = theta[1]
    tau = theta[2]
    r = theta[3]
    b = theta[4]
    if (log==TRUE){
      return (dunif(sigma,0,10,log) + dunif(tau,0,10,log) + dunif(r,-10,10,log) + dunif(b,0,10,log))
    }
    else{
      return (dunif(sigma,0,10,log) * dunif(tau,0,10,log) * dunif(r,-10,10,log) * dunif(b,0,10,log))
    }
  }
  # sampler from the initial distribution of the states
  model.kangarooLogistic$rinitial = function(theta,N){
    return (matrix(rlnorm(N,meanlog = 5,sdlog = sqrt(10)),nrow = N))
  }
  # sampler from the transition distribution of the states
  model.kangarooLogistic$rtransition = function(Xt,t,theta){
    sigma = theta[1]
    r = theta[3]
    b = theta[4]
    N = nrow(Xt)
    logXtold = log(Xt)
    dt = data_kangaroo["time",t] - data_kangaroo["time",t-1]
    M = dt/0.001
    delta = dt/M #Euler method discretization step size
    for (i in 1:M) {
      logXtnew = logXtold + (r-b*exp(logXtold))*delta + sigma*sqrt(delta)*rnorm(N)
      logXtold = logXtnew
    }
    return (matrix(exp(logXtnew),nrow = N))
  }
  # density of the observations
  model.kangarooLogistic$dobs = function(Yt,Xt,t,theta,log = TRUE){
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
  model.kangarooLogistic$robs = function(Xt,t,theta){
    tau = theta[2]
    n = 1/tau
    return (rnbinom(2,size = n,mu = Xt))
  }
  # OPTIONAL: lower and upper bounds of observations
  model.kangarooLogistic$lower = c(0,0)
  model.kangarooLogistic$upper = c(Inf,Inf)
  return (model.kangarooLogistic)
}
