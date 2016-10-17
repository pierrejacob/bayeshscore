#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# This file contains some basic Kalman filter (no smoother here)
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

#this function runs a basic kalman filter and returns the predictive means and variances
#of both states and observations, and the filtering means and variances
#Prototype: theta = (phi, psi, sigmaV2, sigmaW2)
kalmanfilter = function(Y,lineargaussianmodel) {
  nobserv = ncol(Y)
  phi = lineargaussianmodel$theta[1]
  psi = lineargaussianmodel$theta[2]
  sigmaV2 = lineargaussianmodel$theta[3]
  sigmaW2 = lineargaussianmodel$theta[4]
  #initialize containers
  muX_t_t_1 = vector("numeric",nobserv) #contains the means of X_t given Y_1,...,Y_t-1
  muX_t_t = vector("numeric",nobserv) #contains the means of X_t given Y_1,...,Y_t
  PX_t_t_1 = vector("numeric",nobserv) #contains the variances of X_t given Y_1,...,Y_t-1
  PX_t_t = vector("numeric",nobserv) #contains the variances of X_t given Y_1,...,Y_t
  muY_t_t_1 = vector("numeric",nobserv) #contains the means of Y_t given Y_1,...,Y_t-1
  PY_t_t_1 = vector("numeric",nobserv) #contains the variances of Y_t given Y_1,...,Y_t-1
  #initialize recursion
  muX_t_t_1[1] = lineargaussianmodel$initialmean
  PX_t_t_1[1] = lineargaussianmodel$initialvar
  Kt = PX_t_t_1[1]*psi/(psi*PX_t_t_1[1]*psi + sigmaV2)
  muX_t_t[1] = muX_t_t_1[1] + Kt*(Y[,1] - psi*muX_t_t_1[1])
  PX_t_t[1] = (1-Kt*psi)*PX_t_t_1[1]
  muY_t_t_1[1] = psi*muX_t_t_1[1]
  PY_t_t_1[1] = (psi^2)*PX_t_t_1[1] + sigmaV2
  #iterate
  for (t in 2:nobserv) {
    muX_t_t_1[t] = phi*muX_t_t[t-1]
    PX_t_t_1[t] = (phi^2)*PX_t_t[t-1] + sigmaW2
    Kt = PX_t_t_1[t]*psi/(psi*PX_t_t_1[t]*psi + sigmaV2)
    muX_t_t[t] = muX_t_t_1[t] + Kt*(Y[,t] - psi*muX_t_t_1[t])
    PX_t_t[t] = (1-Kt*psi)*PX_t_t_1[t]
    muY_t_t_1[t] = psi*muX_t_t_1[t]
    PY_t_t_1[t] = (psi^2)*PX_t_t_1[t] + sigmaV2
  }
  return (list(muX_t_t_1 = muX_t_t_1, muX_t_t = muX_t_t, muY_t_t_1 = muY_t_t_1,
               PX_t_t_1 = PX_t_t_1, PX_t_t = PX_t_t, PY_t_t_1 = PY_t_t_1))
}

