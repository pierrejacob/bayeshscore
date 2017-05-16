#'@rdname get_model_poissonhmm
#'@title get_model_poissonhmm
#'@description Poisson Hidden-Markov Model
#'@export
get_model_poissonhmm <- function(){
  model = list()
  model$observation_type = "discrete"
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

  # density of the observations
  model$dobs = function(Yt,Xt,t,theta,log = TRUE){
    l = model$lambda
    lx = model$lambdaX
    return (dpois(Yt,l+lx*Xt,log = log))
  }


  # OPTIONAL: likelihood of the observations from time 1 to t
  # This relies on some Kalman filter (passed as a byproduct)
  model$likelihood = function(observations,t,theta,log = TRUE){
    y_1_t = observations[,1:t]
    if(sum(y_1_t<0)>0){
      if (log) {return(-Inf)}
      else {return(0)}
    } else {
      l = model$lambda
      lx = model$lambdaX
      all_paths = expand.grid(rep(list(0:1), t))
      if (t==1){
        temp = 0
        for (i in 1:(2^t)){
          X = all_paths[i,1:t]
          rate = l+X*lx
          temp = temp + prod(exp(-(rate))*(rate^y_1_t)/factorial(y_1_t))*(1/2) + 0*theta
          #the "0*theta" term (artificially) allows integrate to recognize this as a function of theta
        }
        if (log) {return(log(temp))}
        else {return(temp)}

      } else {
        temp = 0
        for (i in 1:(2^t)){
          X = all_paths[i,1:t]
          rate = l+X*lx
          temp = temp + prod(exp(-(rate))*(rate^y_1_t)/factorial(y_1_t))*(1/2)*(theta^sum((1-X[2:t])*(1-X[1:(t-1)])+X[2:t]*X[1:(t-1)]))*((1-theta)^sum(X[2:t]*(1-X[1:(t-1)])+(1-X[2:t])*X[1:(t-1)]))
        }
        if (log) {return(log(temp))}
        else {return(temp)}
      }
    }
  }

  model$dpredictive = function(observations,t,theta,log = TRUE){
    if (t==1){
      return(model$likelihood(observations,t,theta,log))
    } else {
      if (log){return(model$likelihood(observations,t,theta,log)-model$likelihood(observations,t-1,theta,log))}
      else {return(model$likelihood(observations,t,theta,log)/model$likelihood(observations,t-1,theta,log))}
    }
  }

  # OPTIONAL: simulate observations
  model$robs = function(Xt,t,theta){
    l = model$lambda
    lx = model$lambdaX
    N = ncol(Xt)
    return (matrix(rpois(N,l+lx*Xt),ncol=N))
  }

  # OPTIONAL: lower and upper bounds of observations
  model$lower = c(0)
  model$upper = c(Inf)
  return(model)
}
