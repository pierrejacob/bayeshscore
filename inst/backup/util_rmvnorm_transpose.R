#'@rdname fast_rmvnorm_transpose
#'@title fast_rmvnorm_transpose
#'@description fast_rmvnorm_transpose
#'@export

fast_rmvnorm_transpose <- function(nparticles, mean, covariance){
  return(rmvnorm_transpose(nparticles, mean, covariance))
}

#'@rdname fast_rmvnorm_transpose_cholesky
#'@title fast_rmvnorm_transpose_cholesky
#'@description fast_rmvnorm_transpose_cholesky
#'@export

fast_rmvnorm_transpose_cholesky <- function(nparticles, mean, cholcovariance){
  return(rmvnorm_transpose_cholesky(nparticles, mean, cholcovariance))
}
