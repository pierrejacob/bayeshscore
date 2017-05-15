#'@rdname get_model_poissonhmm
#'@title get_model_poissonhmm
#'@description Poisson Hidden-Markov Model
#'@export
get_model_poissonhmm <- function(){
  model = list()
  model$lambda = 1
  model$lambdaX = 3

  # dimension of parameter
  model$dimtheta = 1
  model$dimY = 1
  model$dimX = 1

  # sampler from the prior distribution on parameters
  model$rprior = function(Ntheta){
    return (rbind(runif(Ntheta,0,1)))
  }

  # prior distribution density on parameters
  model$dprior = function(theta, log = TRUE){
    return (dunif(theta,min = 0,max = 1,log = log))
  }

  # sampler from the initial distribution of the states
  model$rinitial = function(theta,N){
    return (matrix(sample(c(0,1),N,replace = TRUE), ncol = N))
  }

  # sampler from the transition density of the states
  model$rtransition = function(Xt,t,theta){
    Xnew = Xt
    N = ncol(Xnew)
    index_0 = (Xnew == 0)
    N0 = sum(index_0)
    Xnew[,index_0] = sample(c(0,1),N0,replace = TRUE,prob = c(theta,1-theta))
    Xnew[,!index_0] = sample(c(0,1),N-N0,replace = TRUE,prob = c(1-theta,theta))
    return (matrix(Xnew, ncol = N))
  }


  ############################################################################################


  # density of the observations
  model$dobs = function(Yt,Xt,t,theta,log = TRUE){
    l = model$lambda
    lx = model$lambdaX
    return (dpois(Yt,l+lx*Xt,log = log))
  }

  # OPTIONAL: simulate observations
  model$robs = function(Xt,t,theta){
    l = model$lambda
    lx = model$lambdaX
    N = nrow(Xt)
    return (rpois(N,l+lx*Xt))
  }

  # OPTIONAL: lower and upper bounds of observations
  model$lower = c(0)
  model$upper = c(Inf)
  return(model)
}
