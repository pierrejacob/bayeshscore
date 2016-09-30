#'@export
rtransition_logistic <- function(Xt, t, theta){
  return(matrix(rtransition_logistic_c(Xt[,1], data_kangaroo["time",t] - data_kangaroo["time",t-1], theta[1], theta[3], theta[4]), ncol = 1))
}
