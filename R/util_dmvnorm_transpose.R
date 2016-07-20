#'@rdname fast_dmvnorm_transpose
#'@title fast_dmvnorm_transpose
#'@description fast_dmvnorm_transpose
#'@export

fast_dmvnorm_transpose <- function(x, mean, chol_covariance){
  return(dmvnorm_transpose(x, mean, chol_covariance))
}

#'@rdname fast_dmvnorm_transpose_cholesky
#'@title fast_dmvnorm_transpose_cholesky
#'@description fast_dmvnorm_transpose_cholesky
#'@export
fast_dmvnorm_transpose_cholesky <- function(x, mean, cholcovariance){
  return(dmvnorm_transpose_cholesky(x, mean, cholcovariance))
}
