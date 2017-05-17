#'@rdname get_batch_initial_proposal
#'@title get_batch_initial_proposal
#'@description This function provides an initial proposal targetting P(theta|Ybatch) based on fitting a Normal to a PMMH sample from P(theta|Ybatch)
#'@export
get_batch_initial_proposal = function(observations,model,algorithmic_parameters){
  dimtheta = model$dimtheta
  #extract batch of observations
  b = algorithmic_parameters$initialbatchsize
  observations_1_to_b = matrix(observations[,1:b], nrow = model$dimY)
  #approximate posterior P(theta|batch) by a Normal
  PMCMCsample = PMMH(observations_1_to_b,model,algorithmic_parameters)
  covariance = cov.wt(PMCMCsample)
  mean_b <- covariance$center
  cov_b <- matrix(covariance$cov,nrow = dimtheta) + diag(rep(10^(-4)/dimtheta),dimtheta)
  #define corresponding proposal for theta
  rinitial = function(Ntheta) {
    return (fast_rmvnorm(Ntheta,mean_b,cov_b))
  }
  logdinitial = function(theta, log = TRUE) {
    if (is.null(dim(theta))) {
      if (log) {
        return (fast_dmvnorm(matrix(theta,ncol = model$dimtheta),mean_b,cov_b))
      }
      else {
        return (exp(fast_dmvnorm(matrix(theta,ncol = model$dimtheta),mean_b,cov_b)))
      }
    }
    else {
      if (log) {
        return (fast_dmvnorm(theta,mean_b,cov_b))
      }
      else {
        return (exp(fast_dmvnorm(theta,mean_b,cov_b)))
      }
    }
  }
  return (list(mean = mean_b, cov = cov_b, r = rinitial, logd = logdinitial))
}
