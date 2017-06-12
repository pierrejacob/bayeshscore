#'@rdname get_model_lineargaussiansimpler_discreteprior
#'@title get_model_lineargaussiansimpler_discreteprior
#'@description Univariate linear Gaussian model with one parameter (\code{sigmaV2})and discrete prior
#'Latent states: X[t] = \code{phi}*X[t-1] + N(0,\code{sigmaW2}).
#'Observations: Y[t] = \code{psi}*X[t] + N(0,\code{sigmaV2}).
#'@export
get_model_lineargaussiansimpler_discreteprior <- function(){
  model = list()
  model$observation_type = 'continuous'
  model$supportprior = seq(0.1,3,0.01)
  model$phi = 0.5
  model$psi = 0.5
  model$sigmaW2 = 1

  # dimension of parameter
  model$dimtheta = 1
  model$dimY = 1
  model$dimX = 1

  # sampler from the prior distribution on parameters
  model$rprior = function(Ntheta){
    support = model$supportprior
    sigmaV2 = sample(support,Ntheta,replace = TRUE)
    return (rbind(sigmaV2))
  }

  # prior distribution density on parameters
  model$dprior = function(theta, log = TRUE){
    if (theta %in% model$supportprior) {
      n = length(model$supportprior)
      if (log) {return (-log(n))}
      else {return (1/n)}
    }
    else {
      if (log) {return (-Inf)}
      else {return (0)}
    }
  }

  # sampler from the initial distribution of the states
  model$rinitial = function(theta,N){
    phi = model$phi
    sigmaW2 = model$sigmaW2
    return (matrix(rnorm(N, mean = 0, sd = sqrt((sigmaW2)/(1-phi^2))), ncol = N))
  }

  # sampler from the transition density of the states
  model$rtransition = function(Xs,t,theta){
    phi = model$phi
    sigmaW2 = model$sigmaW2
    N = ncol(Xs)
    return (matrix(phi*Xs + rnorm(N, mean = 0, sd = sqrt(sigmaW2)), ncol = N))
  }

  # density of the observations
  model$dobs = function(Yt,Xt,t,theta,log = TRUE){
    psi = model$psi
    sigmaV2 = theta[1]
    return (dnorm(Yt,mean = psi*Xt,sd = sqrt(sigmaV2), log))
  }

  # first and second partial derivatives of the observation log-density
  # The function is vectorized with respect to the states Xt (dimX by Nx), so that it outputs:
  # >> the jacobian (Nx by dimY matrix: each row is the transpose of the corresponding gradients row-wise)
  # >> the Hessian diagonals (Nx by dimY matrix: each row is the diagonal coeffs of the corresponding Hessian)
  model$derivativelogdobs = function(Yt,Xt,t,theta){
    psi = model$psi
    sigmaV2 = theta[1]
    N = ncol(Xt)
    d1 = t((psi*Xt-repeat_column(N,Yt))/sigmaV2)
    d2 = matrix(-1/sigmaV2,nrow = N, ncol = model$dimY)
    return (list(jacobian = d1, hessiandiag = d2))
  }

  # OPTIONAL: likelihood of the observations from time 1 to t
  # This relies on some Kalman filter (passed as a byproduct)
  model$likelihood = function(observations,t,theta,KF,log = TRUE){
    phi = model$phi
    sigmaW2 = model$sigmaW2
    psi = model$psi
    sigmaV2 = theta[1]
    initial_mean = 0
    initial_var = (sigmaW2)/(1-phi^2)
    # we make the likelihood an explicit function of the observation at time t
    # to allow the computation of the derivative of the log-predictive density
    KF = KF_assimilate_one(observations[,t,drop=FALSE],t,phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var,KF)
    ll = 0
    for (i in 1:t){
      ll = ll + KF_logdpredictive(observations[,i,drop=FALSE],i,KF)
    }
    if (log) {return(ll)}
    else {return(exp(ll))}
  }

  # OPTIONAL: one-step predicitve density of the observation at time t given all the past from 1 to (t-1)
  # This relies on some Kalman filter (passed as a byproduct)
  model$dpredictive = function(observations,t,theta,KF,log = TRUE){
    phi = model$phi
    sigmaW2 = model$sigmaW2
    psi = model$psi
    sigmaV2 = theta[1]
    initial_mean = 0
    initial_var = (sigmaW2)/(1-phi^2)
    # we make it an explicit function of the observation at time t to allow computation of the derivative
    KF = KF_assimilate_one(observations[,t,drop=FALSE],t,phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var,KF)
    incremental_ll = KF_logdpredictive(observations[,t,drop=FALSE],t,KF)
    if (log) {return(incremental_ll)}
    else {return(exp(incremental_ll))}
  }

  # OPTIONAL: derivatives of the predicitve density
  # The function outputs:
  # >> the transpose of the gradient (1 by dimY)
  # >> the Hessian diagonal coefficients (1 by dimY)
  model$derivativelogdpredictive = function(observations,t,theta,KF,dimY) {
    m = KF[[t]]$muY_t_t_1
    V = KF[[t]]$PY_t_t_1
    d1 = t(matrix((m-observations[,t])/V))
    d2 = matrix(-1/V)
    return (list(jacobian = d1, hessiandiag = d2))
  }

  # OPTIONAL: initialize byproducts (e.g. Kalman filters, etc ...)
  model$initialize_byproducts = function(theta, observations){
    KF = vector("list",ncol(observations))
    return(KF)
  }

  # OPTIONAL: update byproducts (e.g. Kalman filters, etc ...)
  model$update_byproduct = function(KF, t, theta, observations){
    phi = model$phi
    sigmaW2 = model$sigmaW2
    psi = model$psi
    sigmaV2 = theta[1]
    initial_mean = 0
    initial_var = (sigmaW2)/(1-phi^2)
    KF_updated = KF_assimilate_one(observations[,t,drop=FALSE],t,phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var, KF)
    return(KF_updated)
  }

  # OPTIONAL: simulate observations
  model$robs = function(Xt,t,theta){
    psi = model$psi
    sigmaV2 = theta[1]
    N = ncol(Xt)
    return (matrix(psi*Xt + rnorm(N, mean = 0, sd = sqrt(sigmaV2)),ncol = N))
  }

  return(model)
}
