#'@export
rtransition_logistic <- function(Xt, t, theta,timesteps){
  return(matrix(rtransition_logistic_c(Xt[,1], timesteps[t] - timesteps[t-1], theta[1], theta[3], theta[4]), ncol = 1))
}
