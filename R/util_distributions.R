#------------------------------------------------------------------------------#
#-----------------  Some useful distributions (density and sampler) -----------#
#------------------------------------------------------------------------------#
#'@rdname fast_dmvnorm
#'@title fast_dmvnorm
#'@description Fast Normal log-density using Rcpp (observations row-wise, dimensions column-wise)
#'@export
fast_dmvnorm <- function(x, mean, covariance){
  return(dmvnorm(x, mean, covariance))
}
#'@rdname fast_rmvnorm
#'@title fast_rmvnorm
#'@description Fast samples from Normal using Rcpp (observations row-wise, dimensions column-wise)
#'@export
fast_rmvnorm <- function(nparticles, mean, covariance){
  return(rmvnorm(nparticles, mean, covariance))
}
#------------------------------------------------------------------------------#
#'@rdname fast_dmvnorm_transpose
#'@title fast_dmvnorm_transpose
#'@description Fast Normal log-density using Rcpp (dimensions row-wise, observations column-wise)
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
#'@rdname fast_rmvnorm_transpose
#'@title fast_rmvnorm_transpose
#'@description Fast samples from  Normal using Rcpp (dimensions row-wise, observations column-wise)
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
#------------------------------------------------------------------------------#
#'@rdname dinvgamma
#'@title dinvgamma
#'@description Density of inverse gamma distribution
#'@export
dinvgamma = function(y,a,b,log){
  is_positive = (y>0)
  result = y
  y_pos = y[is_positive]
  y_neg = y[!is_positive]
  result[is_positive] = (a*log(b)-lgamma(a)-(a+1)*log(y_pos)-b/y_pos)
  result[!is_positive] = -Inf
  if (log==TRUE) {
    return (result)
  }
  else {
    return(exp(result))
  }
}
# # Sanity check: density of inverse Gamma distribution
# a = 10
# b = 1
# g11 <- ggplot(data.frame(sim=rinvgamma(10000,a,b)), aes(x = sim)) + geom_density() +
#   stat_function(fun = function(y)dinvgamma(y,a,b,FALSE), colour = "red",size=1.5,linetype=2)
# grid.arrange(g11,ncol = 1, nrow = 1)

#'@rdname rinvgamma
#'@title rinvgamma
#'@description Samples from inverse gamma distribution
#'@export
rinvgamma = function(N,a,b){
  return (1/rgamma(N,a,rate = b))
}
# # Sanity check: density of inverse Gamma distribution
# a = 10
# b = 1
# g11 <- ggplot(data.frame(sim=rinvgamma(10000,a,b)), aes(x = sim)) + geom_density() +
#   stat_function(fun = function(y)dinvgamma(y,a,b,FALSE), colour = "red",size=1.5,linetype=2)
# grid.arrange(g11,ncol = 1, nrow = 1)
#------------------------------------------------------------------------------#
#'@rdname dinvchisq
#'@title dinvchisq
#'@description Density of scaled inverse chi-square distribution
#'@export
dinvchisq = function(y,df,s2=1,log){
  return (dinvgamma(y,a = df/2,b = df*s2/2,log))
}
# # Sanity check: density of scaled inverse chi-square distribution
# df = 20
# s2 = 2
# g11 <- ggplot(data.frame(sim=rinvchisq(10000,df,s2)), aes(x = sim)) + geom_density() +
#   stat_function(fun = function(y)dinvchisq(y,df,s2,FALSE), colour = "red",size=1.5,linetype=2)
# grid.arrange(g11,ncol = 1, nrow = 1)

#'@rdname rinvchisq
#'@title rinvchisq
#'@description Samples from scaled inverse chi-square distribution
#'@export
rinvchisq = function(N,df,s2=1){
  return (rinvgamma(N,a = df/2,b = df*s2/2))
}
# # Sanity check: density of scaled inverse chi-square distribution
# df = 20
# s2 = 2
# g11 <- ggplot(data.frame(sim=rinvchisq(10000,df,s2)), aes(x = sim)) + geom_density() +
#   stat_function(fun = function(y)dinvchisq(y,df,s2,FALSE), colour = "red",size=1.5,linetype=2)
# grid.arrange(g11,ncol = 1, nrow = 1)
#------------------------------------------------------------------------------#
#'@rdname dtscaled
#'@title dtscaled
#'@description Density of scaled Student t distribution
#'@export
dtscaled = function(y,df,s2,log){
  if (log==TRUE) {
    return (lgamma((df+1)/2)-lgamma(df/2)-(1/2)*log(df*pi*s2)-((df+1)/2)*log(1+(1/df)*(y)^2/s2))
  }
  else {
    return (exp(lgamma((df+1)/2)-lgamma(df/2)-(1/2)*log(df*pi*s2)-((df+1)/2)*log(1+(1/df)*(y)^2/s2)))
  }
}
# Sanity check: density of scaled t distribution
# df = 5
# s2 = 2
# g11 <- ggplot(data.frame(sim=rtscaled(10000,df,s2)), aes(x = sim)) + geom_density() +
#   stat_function(fun = function(y)dtscaled(y,df,s2,FALSE), colour = "red",size=1.5,linetype=2)
# grid.arrange(g11,ncol = 1, nrow = 1)

#'@rdname rtscaled
#'@title rtscaled
#'@description Samples from scaled Student t distribution
#'@export
rtscaled = function(N,df,s2){
  return (sqrt(s2)*rt(N,df))
}
# Sanity check: density of scaled t distribution
# df = 5
# s2 = 2
# g11 <- ggplot(data.frame(sim=rtscaled(10000,df,s2)), aes(x = sim)) + geom_density() +
#   stat_function(fun = function(y)dtscaled(y,df,s2,FALSE), colour = "red",size=1.5,linetype=2)
# grid.arrange(g11,ncol = 1, nrow = 1)
#------------------------------------------------------------------------------#


