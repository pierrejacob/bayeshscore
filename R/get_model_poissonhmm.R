#'@rdname get_model_poissonhmm
#'@title get_model_poissonhmm
#'@description This implements a Poisson HMM
#'@export
get_model_poissonhmm <- function(){
  model.poissonhmm = list()
  model.poissonhmm$lambda = 1
  model.poissonhmm$lambdaX = 3

  # dimension of parameter
  model.poissonhmm$dimtheta = 1
  model.poissonhmm$dimY = 1
  model.poissonhmm$dimX = 1

  # sampler from the prior distribution on parameters
  model.poissonhmm$rprior = function(Ntheta){
    return (cbind(runif(Ntheta,0,1)))
  }

  # prior distribution density on parameters
  model.poissonhmm$dprior = function(theta, log = TRUE){
    return (dunif(theta,min = 0,max = 1,log = log))
  }

  # sampler from the initial distribution of the states
  model.poissonhmm$rinitial = function(theta,N){
    return (matrix(sample(c(0,1),N,replace = TRUE), nrow = N))
  }

  # sampler from the transition density of the states
  model.poissonhmm$rtransition = function(Xt,t,theta){
    Xnew = Xt
    N = nrow(Xnew)
    index_0 = (Xnew == 0)
    N0 = sum(index_0)
    Xnew[index_0,] = sample(c(0,1),N0,replace = TRUE,prob = c(theta,1-theta))
    Xnew[!index_0,] = sample(c(0,1),N-N0,replace = TRUE,prob = c(1-theta,theta))
    return (matrix(Xnew, nrow = N))
  }


  ############################################################################################


  # density of the observations
  model.poissonhmm$dobs = function(Yt,Xt,t,theta,log = TRUE){
    l = model.poissonhmm$lambda
    lx = model.poissonhmm$lambdaX
    return (dpois(Yt,l+lx*Xt,log = log))
  }

  # OPTIONAL: simulate observations
  model.poissonhmm$robs = function(Xt,t,theta){
    l = model.poissonhmm$lambda
    lx = model.poissonhmm$lambdaX
    N = nrow(Xt)
    return (rpois(N,l+lx*Xt))
  }

  # OPTIONAL: lower and upper bounds of observations
  model.poissonhmm$lower = c(0)
  model.poissonhmm$upper = c(Inf)
  return(model.poissonhmm)
}
