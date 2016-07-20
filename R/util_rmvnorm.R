#'@rdname fast_rmvnorm
#'@title fast_rmvnorm
#'@description fast_rmvnorm
#'@export

fast_rmvnorm <- function(nparticles, mean, covariance){
  return(rmvnorm(nparticles, mean, covariance))
}
