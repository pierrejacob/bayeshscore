#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------- Some functions for Kalman filter (no smoother here) -----------#
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#'@rdname KF_assimilate_one
#'@title KF_assimilate_one
#'@description This function performs one step of Kalman filtering
#'@export
KF_assimilate_one = function(Yt, t, phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var, KF_current) {
  # The Rcpp implementation takes matrices as arguments
  if (is.null(dim(phi))) {phi = matrix(phi)}
  if (is.null(dim(psi))) {psi = matrix(psi)}
  if (is.null(dim(sigmaV2))) {sigmaV2 = matrix(sigmaV2)}
  if (is.null(dim(sigmaW2))) {sigmaW2 = matrix(sigmaW2)}
  if (is.null(dim(initial_mean))) {initial_mean = matrix(initial_mean)}
  if (is.null(dim(initial_var))) {initial_var = matrix(initial_var)}
  # Note: index in Cpp starts at 0, hence the time shift (t-1)
  return (KF_assimilate_one_cpp(Yt, t-1,initial_mean, initial_var, phi, psi, sigmaV2, sigmaW2, KF_current))
}
#----------------------------------------------------------------------------------------
#'@rdname KF_filtering
#'@title KF_filtering
#'@description This function runs a basic kalman filter and returns the predictive means and variances
#'of both states and observations, and the filtering means and variances
#'Prototype: state parameters = phi, sigmaV2 // observation parameter: psi, sigmaW2
#'@export
KF_filtering = function(Y,phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var) {
  # The Rcpp implementation takes matrices as arguments
  if (is.null(dim(phi))) {phi = matrix(phi)}
  if (is.null(dim(psi))) {psi = matrix(psi)}
  if (is.null(dim(sigmaV2))) {sigmaV2 = matrix(sigmaV2)}
  if (is.null(dim(sigmaW2))) {sigmaW2 = matrix(sigmaW2)}
  if (is.null(dim(initial_mean))) {initial_mean = matrix(initial_mean)}
  if (is.null(dim(initial_var))) {initial_var = matrix(initial_var)}
  return (KF_filtering_cpp(Y, initial_mean, initial_var, phi, psi, sigmaV2, sigmaW2))
}
#----------------------------------------------------------------------------------------
#'@rdname KF_logdpredictive
#'@title KF_logdpredictive
#'@description This function returns the predictive log-density of the next observation
#'@export
KF_logdpredictive = function(Yt, t, KF_current) {
  m = matrix(KF_current[[t]]$muY_t_t_1, nrow = nrow(Yt), ncol = 1)
  V = matrix(KF_current[[t]]$PY_t_t_1,nrow = nrow(Yt), ncol = nrow(Yt))
  return(fast_dmvnorm_transpose(Yt,m,V))
}

#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# Below is the initial implementation in R. The final implementation was done in Rcpp
# to improve computation speed.
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# KF_assimilate_one_R= function(Yt, t, phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var, KF_current) {
#   if (t==1) {
#     muX_t_t_1 = initial_mean
#     PX_t_t_1 = initial_var
#     Kt = PX_t_t_1%*%t(psi)%*%solve(psi%*%PX_t_t_1%*%t(psi) + sigmaV2)
#     muX_t_t = muX_t_t_1 + Kt%*%(Yt - psi%*%muX_t_t_1)
#     if (!is.null(dim(PX_t_t_1))) {
#       PX_t_t = (diag(1,nrow(PX_t_t_1))-Kt%*%psi)%*%PX_t_t_1
#     } else {
#       PX_t_t = (1-Kt%*%psi)%*%PX_t_t_1
#     }
#     muY_t_t_1 = psi%*%muX_t_t_1
#     PY_t_t_1 = psi%*%PX_t_t_1%*%t(psi) + sigmaV2
#   } else {
#     muX_t_t_1 = phi%*%KF_current[[t-1]]$muX_t_t
#     PX_t_t_1 = phi%*%KF_current[[t-1]]$PX_t_t%*%t(phi) + sigmaW2
#     Kt = PX_t_t_1%*%t(psi)%*%solve(psi%*%PX_t_t_1%*%t(psi) + sigmaV2)
#     muX_t_t = muX_t_t_1 + Kt%*%(Yt - psi%*%muX_t_t_1)
#     if (!is.null(dim(PX_t_t_1))) {
#       PX_t_t = (diag(1,nrow(PX_t_t_1))-Kt%*%psi)%*%PX_t_t_1
#     } else {
#       PX_t_t = (1-Kt%*%psi)%*%PX_t_t_1
#     }
#     muY_t_t_1 = psi%*%muX_t_t_1
#     PY_t_t_1 = psi%*%PX_t_t_1%*%t(psi) + sigmaV2
#   }
#   KF_updated = KF_current
#   KF_updated[[t]] = list(muX_t_t_1 = muX_t_t_1, #contains the means of X_t given Y_1,...,Y_t-1
#                          PX_t_t_1 = PX_t_t_1, #contains the variances of X_t given Y_1,...,Y_t-1
#                          muX_t_t = muX_t_t, #contains the means of X_t given Y_1,...,Y_t
#                          PX_t_t = PX_t_t, #contains the variances of X_t given Y_1,...,Y_t
#                          muY_t_t_1 = muY_t_t_1, #contains the means of Y_t given Y_1,...,Y_t-1
#                          PY_t_t_1 = PY_t_t_1) #contains the variances of Y_t given Y_1,...,Y_t-1
#   return (KF_updated)
# }
#----------------------------------------------------------------------------------------
# KF_filtering_R = function(Y,phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var) {
#   nobs = ncol(Y)
#   KF = vector("list",nobs)
#   for (t in 1:nobs) {
#     KF = KF_assimilate_one_R(Y[,t,drop=FALSE],t,phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var, KF)
#   }
#   return (KF)
# }
