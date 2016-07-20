#'@rdname fast_dmvnorm
#'@title fast_dmvnorm
#'@description fast_dmvnorm
#'@export

fast_dmvnorm <- function(x, mean, covariance){
  return(dmvnorm(x, mean, covariance))
}
