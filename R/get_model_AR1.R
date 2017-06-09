#'@rdname get_model_AR1
#'@title get_model_AR1
#'@description Auto-regressive model of order 1
#'@export
get_model_AR1 <- function(nu0, sigma02){
  model = list()
  model$observation_type = 'continuous'

  # dimension of parameter
  model$dimtheta = 2
  model$dimY = 1

  # sampler from the prior distribution on parameters
  model$rprior = function(Ntheta){
    phi = runif(Ntheta,-1,1)
    sigma2 = rinvchisq(Ntheta,nu0,sigma02)
    return (rbind(phi,sigma2))
  }

  # prior distribution density on parameters
  model$dprior = function(theta, log = TRUE){
    logd = dunif(theta[1],-1,1,log = TRUE) + dinvchisq(theta[2],nu0,sigma02,log = TRUE)
    if (log==TRUE) {return (logd)}
    else {return (exp(logd))}
  }

  # one-step predicitve density of the observation at time t given all the past from 1 to (t-1)
  model$dpredictive = function(observations,t,theta,log = TRUE){
    if (t==1) {
      mu = 0
      sigma2 = theta[2]/(1-theta[1]^2)
    } else {
      mu = theta[1]*observations[,t-1]
      sigma2 = theta[2]
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
    if (t==1) {
      mu = 0
      sigma2 = theta[2]/(1-theta[1]^2)
    } else {
      mu = theta[1]*observations[,t-1]
      sigma2 = theta[2]
    }
    deriv1 <- -(observations[,t]-mu)/sigma2
    deriv2 <- -1/sigma2
    return (list(jacobian = matrix(deriv1, 1, 1), hessiandiag = matrix(deriv2, 1, 1)))
  }

  # OPTIONAL: simulate observations
  model$robs = function(nobservations,theta){
    Y = matrix(NA, nrow = 1, ncol = nobservations)
    Y[,1] = rnorm(1, mean = 0, sd = sqrt(theta[2]/(1-theta[1]^2)))
    if (nobservations > 1){
      for (t in 2:nobservations){
        Y[,t] = theta[1]*Y[,t-1] + rnorm(1, mean = 0, sd = sqrt(theta[2]))
      }
    }
    return (Y)
  }

  # OPTIONAL: simulate Ny draws of y_t given theta and the past y_1 to y_(t-1)
  # (with the convention y_0 = NULL)
  # outputs: matrix of Ny draws of Yt given theta and past (dimY by Ny matrix)
  model$rpredictive = function(Ny,t,theta,y_past){
    if (t == 1) {
      Yt = matrix(rnorm(Ny, mean = 0, sd = sqrt(theta[2]/(1-theta[1]^2))), ncol = Ny)
    } else if (t >= 2) {
      Yt = theta[1]*repeat_column(Ny,y_past[,t-1,drop=FALSE]) + matrix(rnorm(Ny*model$dimY, mean = 0, sd = sqrt(theta[2])),ncol=Ny)
    }
    return (Yt)
  }


  return(model)
}
