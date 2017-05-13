#'@rdname get_model_simplerlineargaussian
#'@title get_model_simplerlineargaussian
#'@description This implements the univariate linear Gaussian model
#'@export
get_model_simplerlineargaussian <- function(){
  model = list()
  model$psi = 1
  model$sigmaV2 = 0.8

  # dimension of parameter
  model$dimtheta = 2
  model$dimY = 1
  model$dimX = 1

  # sampler from the prior distribution on parameters
  model$rprior = function(Ntheta){
    phi = runif(Ntheta,0.1,0.9)
    sigmaW2 = runif(Ntheta,0.1,10)
    return (rbind(phi,sigmaW2))
  }

  # prior distribution density on parameters
  model$dprior = function(theta, log = TRUE){
    phi = theta[1]
    sigmaW2 = theta[2]
    if (log==TRUE){
      return (dunif(phi,0.1,0.9,log)+dunif(sigmaW2,0.1,10,log))
    }
    else{
      return (dunif(phi,0.1,0.9,log)*dunif(sigmaW2,0.1,10,log))
    }
  }

  # sampler from the initial distribution of the states
  model$rinitial = function(theta,N){
    phi = theta[1]
    sigmaW2 = theta[2]
    return (matrix(rnorm(N, mean = 0, sd = sqrt((sigmaW2)/(1-phi^2))), ncol = N))
  }

  # sampler from the transition density of the states
  model$rtransition = function(Xt,t,theta){
    phi = theta[1]
    sigmaW2 = theta[2]
    N = ncol(Xt)
    return (matrix(phi*Xt + rnorm(N, mean = 0, sd = sqrt(sigmaW2)), ncol = N))
  }

  # density of the observations
  model$dobs = function(Yt,Xt,t,theta,log = TRUE){
    return (dnorm(Yt,mean = model$psi*Xt,sd = sqrt(model$sigmaV2), log))
  }

  # first and second partial derivatives (k-th coordinate) of the observation log-density
  model$derivativelogdobs = function(Yt,Xt,t,theta,dimY){
    N = ncol(Xt)
    d1 = matrix((model$psi*Xt-repeat_column(N,Yt))/model$sigmaV2,nrow = N, ncol = dimY)
    d2 = matrix(-1/model$sigmaV2,nrow = N, ncol = dimY)
    return (list(jacobian = d1, hessiandiag = d2))
  }

  # OPTIONAL: likelihood of the observations from time 1 to t
  # This relies on some Kalman filter (passed as a byproduct)
  model$likelihood = function(observations,t,theta,byproduct,log = TRUE){
    incremental_ll <- byproduct$get_incremental_ll()
    if (log) {
      return(sum(incremental_ll[1:t]))
    } else {
      return(exp(sum(incremental_ll[1:t])))
    }
  }

  # # OPTIONAL: one-step predicitve density of the observation at time t given all the past from time 1 to (t-1)
  # model$dpredictive = function(observations,t,theta,byproduct,log = TRUE){
  #   if (t==1) {
  #     return(model$likelihood(observations,t,theta,byproduct,log))
  #   } else {
  #     if (log) {
  #       return(model$likelihood(observations,t,theta,byproduct,log)-model$likelihood(observations,t-1,theta,byproduct,log))
  #     } else {
  #       return(model$likelihood(observations,t,theta,byproduct,log)/model$likelihood(observations,t-1,theta,byproduct,log))
  #     }
  #   }
  # }

  # # first and second partial derivatives (k-th coordinate) of the one-step predictive log-density
  # model$derivativelogdpredictive = function(observations,t,theta,byproduct){
  #   if (t==1){
  #     logpredictive = function(y) {model$dpredictive(y,t,theta,byproduct,log = TRUE)}
  #   } else {
  #     logpredictive = function(y) {model$dpredictive(cbind(observations[1:(t-1)],y),t,theta,byproduct,log = TRUE)}
  #   }
  #   gradient_logdpredictive = grad(logpredictive,observations[,t])
  #   hessian_logdpredictive = hessian(logpredictive,observations[,t])
  #   return (list(grad = gradient_logdpredictive, hessian = hessian_logdpredictive))
  # }

  # OPTIONAL: initialize byproducts (e.g. Kalman filters, etc ...)
  model$initialize_byproducts = function(theta, observations, Ntheta){
    KF <- new(kalman_module$Kalman)
    KF$set_parameters(list(rho = theta[1], sigma = sqrt(theta[2]), eta = 1, tau = sqrt(0.8)))
    KF$set_observations(matrix(observations, ncol = 1))
    KF$first_step()
    return(KF)
  }
  # OPTIONAL: update byproducts (e.g. Kalman filters, etc ...)
  model$update_byproduct = function(KF, t, thetas, observations){
    KF$filtering_step(t-1)
    return(KF)
  }

  # OPTIONAL: simulate observations
  model$robs = function(Xt,t,theta){
    N = ncol(Xt)
    return (matrix(model$psi*Xt + rnorm(N, mean = 0, sd = sqrt(model$sigmaV2)),ncol = N))
  }

  return(model)
}
