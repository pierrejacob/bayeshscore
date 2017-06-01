#'@rdname get_model_ARMA
#'@title get_model_ARMA
#'@description Auto-regressive-Moving-average model of order (p,q)
#'@export
get_model_ARMA <- function(p,q,nu0,sigma02){
  model = list()
  model$observation_type = 'continuous'

  # dimension of parameter
  r = max(p,q+1)
  model$dimtheta = p + q + 1
  model$dimY = 1
  model$dimX = r

  # set some hyperparameters
  model$initial_mean = matrix(0,nrow=model$dimX,ncol=1)
  model$get_initial_var = function(theta){
    sigma2 = theta[model$dimtheta]
    initial_var = diag(sigma2,model$dimX)
    return (initial_var)
  }
  #-------------------------------------------------------------------------------------
  # Get the linear gaussian state-space model representation of the ARMA model
  # This function extends the vector of ARMA coefficients
  # ARMA_coeffs contains the vector of ARMA coefficients (p AR coeffs, then q MA coeffs)
  model$get_extendedcoeffs = function(ARMA_coeffs){
    if ((p>0)&&(q>0)){
      AR_coeffs = c(ARMA_coeffs[1:p],rep(0,r-p))
      MA_coeffs = c(ARMA_coeffs[(p+1):(p+q)],rep(0,r-1-q))
    } else if ((p>0)&&(q==0)) {
      AR_coeffs = c(ARMA_coeffs[1:p],rep(0,r-p))
      MA_coeffs = c(rep(0,r-1-q))
    } else if ((p==0)&&(q>0)) {
      AR_coeffs = c(rep(0,r-p))
      MA_coeffs = c(ARMA_coeffs[(p+1):(p+q)],rep(0,r-1-q))
    }
    extended_coeffs = list(AR = AR_coeffs, MA = MA_coeffs)
    return (extended_coeffs)
  }
  # This function defines the corresponding PHI matrix in the linear gaussian state-space representation
  # ARMA_coeffs contains the vector of ARMA coefficients (p AR coeffs, then q MA coeffs)
  model$get_PHI = function(ARMA_coeffs){
    if (r>=2) {
      extended_coeffs = model$get_extendedcoeffs(ARMA_coeffs)
      PHI = toeplitz(c(0,1,rep(0,r-2)))
      PHI[lower.tri(PHI)] = 0
      PHI[r,] = extended_coeffs$AR[r:1]
      return (PHI)
    } else if (r==1) {
      PHI = matrix(ARMA_coeffs[1])
      return (PHI)
    }
  }
  # This function defines the corresponding PSI matrix in the linear gaussian state-space representation
  # ARMA_coeffs contains the vector of ARMA coefficients (p AR coeffs, then q MA coeffs)
  model$get_PSI = function(ARMA_coeffs){
    if (r>=2){
      extended_coeffs = model$get_extendedcoeffs(ARMA_coeffs)
      PSI = matrix(c(extended_coeffs$MA[(r-1):1],1),nrow = 1)
    } else if (r==1) {
      PSI = matrix(1)
    }
    return (PSI)
  }
  # Define the corresponding state-noise matrix in the linear gaussian state-space representation
  model$get_SIGMAW2 = function(sigma2){
    SIGMAW2 = matrix(0,ncol = r, nrow = r)
    SIGMAW2[r,r] = sigma2
    return(SIGMAW2)
  }
  # Define the corresponding observation-noise matrix in the linear gaussian state-space representation
  model$SIGMAV2 = matrix(0)
  #-------------------------------------------------------------------------------------

  # range for the uniform prior on the ARMA-parameters
  rangeprior = c(-1,1)
  model$rprior = function(Ntheta){
    if (p > 0){phi = matrix(runif(Ntheta*p,rangeprior[1],rangeprior[2]), nrow = p)}
    else {phi = NULL}
    if (q > 0){theta = matrix(runif(Ntheta*q,rangeprior[1],rangeprior[2]), nrow = q)}
    else {theta = NULL}
    sigma2 = rinvchisq(Ntheta,nu0,sigma02)
    return (rbind(phi,theta,sigma2))
  }

  # prior distribution density on parameters
  model$dprior = function(theta, log = TRUE){
    ARMA_coeffs = theta[1:(model$dimtheta-1)]
    sigma2 = theta[model$dimtheta]
    logd = sum(dunif(ARMA_coeffs,rangeprior[1],rangeprior[2],TRUE)) + dinvchisq(sigma2,nu0,sigma02,TRUE)
    if (log==TRUE) {return (logd)}
    else {return (exp(logd))}
  }

  # sampler from the initial distribution of the states
  # (stationary distribution if AR(1), i.e. p = 1 and q = 0)
  model$rinitial = function(theta,N){
    ARMA_coeffs = theta[1:(model$dimtheta-1)]
    sigma2 = theta[model$dimtheta]
    return (matrix(rnorm(N*r, mean = 0, sd = sqrt(model$get_initial_var(theta)[1,1])), ncol = N))
  }


  # sampler from the transition density of the states
  model$rtransition = function(Xt,t,theta){
    ARMA_coeffs = theta[1:(model$dimtheta-1)]
    sigma2 = theta[model$dimtheta]
    PHI = model$get_PHI(ARMA_coeffs)
    N = ncol(Xt)
    Wt = rnorm(N,0,sqrt(sigma2))
    Xnew = PHI%*%Xt
    Xnew[r,] = Xnew[r,] + Wt
    return (Xnew)
  }

  # density of the observations
  model$dobs = function(Yt,Xt,t,theta,log = TRUE){
    ARMA_coeffs = theta[1:(model$dimtheta-1)]
    PSI = model$get_PSI(ARMA_coeffs)
    return (dnorm(Yt,mean = PSI%*%Xt,sd = sqrt(0), log))
  }

  # OPTIONAL: likelihood of the observations from time 1 to t
  # This relies on some Kalman filter (passed as a byproduct)
  model$likelihood = function(observations,t,theta,KF,log = TRUE){
    ARMA_coeffs = theta[1:(model$dimtheta-1)]
    sigma2 = theta[model$dimtheta]
    phi = model$get_PHI(ARMA_coeffs)
    psi = model$get_PSI(ARMA_coeffs)
    sigmaW2 = model$get_SIGMAW2(sigma2)
    sigmaV2 = model$SIGMAV2
    initial_mean = model$initial_mean
    initial_var = model$get_initial_var(theta)
    KF = KF_assimilate_one(observations[,t,drop=FALSE],t,phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var,KF)
    # we make the likelihood an explicit function of the observation at time t
    # to allow the computation of the derivative of the log-predictive density
    ll = 0
    for (i in 1:t){
      ll = ll + KF_logdpredictive(observations[,i,drop=FALSE],i,KF)
    }
    if (log) {
      return(ll)
    } else {
      return(exp(ll))
    }
  }


  # OPTIONAL: one-step predicitve density of the observation at time t given all the past from 1 to (t-1)
  # This relies on some Kalman filter (passed as a byproduct)
  model$dpredictive = function(observations,t,theta,KF,log = TRUE){
    ARMA_coeffs = theta[1:(model$dimtheta-1)]
    sigma2 = theta[model$dimtheta]
    phi = model$get_PHI(ARMA_coeffs)
    psi = model$get_PSI(ARMA_coeffs)
    sigmaW2 = model$get_SIGMAW2(sigma2)
    sigmaV2 = model$SIGMAV2
    initial_mean = model$initial_mean
    initial_var = model$get_initial_var(theta)
    # we make it an explicit function of the observation at time t to allow computation of the derivative
    KF = KF_assimilate_one(observations[,t,drop=FALSE],t,phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var,KF)
    incremental_ll = KF_logdpredictive(observations[,t,drop=FALSE],t,KF)
    if (log) {
      return(incremental_ll)
    } else {
      return(exp(incremental_ll))
    }
  }

  # OPTIONAL: derivatives of the predicitve density
  # The function outputs:
  # >> the transpose of the gradient (1 by dimY)
  # >> the Hessian diagonal coefficients (1 by dimY)
  model$derivativelogdpredictive = function(observations,t,theta,KF) {
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
    ARMA_coeffs = theta[1:(model$dimtheta-1)]
    sigma2 = theta[model$dimtheta]
    phi = model$get_PHI(ARMA_coeffs)
    psi = model$get_PSI(ARMA_coeffs)
    sigmaW2 = model$get_SIGMAW2(sigma2)
    sigmaV2 = model$SIGMAV2
    initial_mean = model$initial_mean
    initial_var = model$get_initial_var(theta)
    KF_updated = KF_assimilate_one(observations[,t,drop=FALSE],t,phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var, KF)
    return(KF_updated)
  }

  # OPTIONAL: simulate observations
  model$robs = function(Xt,t,theta){
    ARMA_coeffs = theta[1:(model$dimtheta-1)]
    PSI = model$get_PSI(ARMA_coeffs)
    return (PSI%*%Xt)
  }

  return(model)
}
