#'@rdname get_model_AR2
#'@title get_model_AR2
#'@description Auto-regressive model of order 2
#'@export
get_model_AR2 <- function(nu0, sigma02){
  model = list()
  model$observation_type = 'continuous'

  # dimension of parameter
  model$dimtheta = 3
  model$dimY = 1


  # uniform distribution on the triangle that guarantees causality of the AR(2) process
  model$runiftriangle = function(){
    phi1 = runif(1,-2,2)
    phi2 = runif(1,-1,1)
    accept = (-(1-phi2) < phi1)&&(phi1 < 1-phi2)
    while (!accept){
      phi1 = runif(1,-2,2)
      phi2 = runif(1,-1,1)
      accept = (-(1-phi2) < phi1)&&(phi1 < 1-phi2)
    }
    return (rbind(phi1,phi2))
  }

  model$rprior = function(Ntheta){
    phi = sapply(1:Ntheta,function(i)model$runiftriangle())
    sigma2 = rinvchisq(Ntheta,nu0,sigma02)
    return (rbind(phi,sigma2))
  }

  # prior distribution density on parameters
  model$dprior = function(theta, log = TRUE){
    if ((-1<theta[2])&&(theta[2]<1)&&(-(1-theta[2])<theta[1])&&(theta[1]<(1-theta[2]))){
      logd = log(1/4) + dinvchisq(theta[3],nu0,sigma02,log = TRUE)
    } else {
      logd = -Inf
    }
    if (log==TRUE) {return (logd)}
    else {return (exp(logd))}
  }

  # one-step predicitve density of the observation at time t given all the past from 1 to (t-1)
  model$dpredictive = function(observations,t,theta,log = TRUE){
    if ((t==1)||(t==2)) {
      mu = 0
      sigma2 = theta[3]/(1-theta[2]^2-(theta[1]^2)*(1+theta[2])/(1-theta[2]))
    } else {
      mu = theta[1]*observations[,t-1]+theta[2]*observations[,t-2]
      sigma2 = theta[3]
    }
    return (dnorm(observations[,t],mean = mu,sd = sqrt(sigma2), log))
  }

  # OPTIONAL: derivatives of the predicitve density with respect to the observation at time t
  # inputs: observations (dimY by T matrix, with T >= t), time index t (int), theta (single vector),
  #         byproduct (OPTIONAL: auxiliary object needed to compute likelihood, e.g. Kalman filter)
  # outputs: list with the following fields
  # jacobian >> the transpose of the gradient (1 by dimY)
  # hessiandiag >> the Hessian diagonal coefficients (1 by dimY)
  # NB: if missing, this field is automatically filled with numerical derivatives
  # via set_default_model in util_default.R)
  model$derivativelogdpredictive = function(observations,t,theta,byproduct) {
    if ((t==1)||(t==2)) {
      mu = 0
      sigma2 = theta[3]/(1-theta[2]^2-(theta[1]^2)*(1+theta[2])/(1-theta[2]))
    } else {
      mu = theta[1]*observations[,t-1]+theta[2]*observations[,t-2]
      sigma2 = theta[3]
    }
    deriv1 <- -(observations[,t]-mu)/sigma2
    deriv2 <- -1/sigma2
    return (list(jacobian = matrix(deriv1, 1, 1), hessiandiag = matrix(deriv2, 1, 1)))
  }

  # OPTIONAL: simulate observations
  model$robs = function(nobservations,theta){
    Y = matrix(NA, nrow = 1, ncol = nobservations)
    Y[,1] = rnorm(1, mean = 0, sd = sqrt(theta[3]/(1-theta[2]^2-(theta[1]^2)*(1+theta[2])/(1-theta[2]))))
    Y[,2] = rnorm(1, mean = 0, sd = sqrt(theta[3]/(1-theta[2]^2-(theta[1]^2)*(1+theta[2])/(1-theta[2]))))
    if (nobservations > 2){
      for (t in 3:nobservations){
        Y[,t] = theta[1]*Y[,t-1] + theta[2]*Y[,t-2] + rnorm(1, mean = 0, sd = sqrt(theta[3]))
      }
    }
    return (Y)
  }
  return(model)
}
