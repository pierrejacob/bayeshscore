#------------------------------------------------------------------------------#
#-------------- Some proposals for the rejuvenation steps of SMC and SMC2  ----#
#------------------------------------------------------------------------------#
#'@rdname get_proposal_independent_normal
#'@title get_proposal_independent_normal
#'@description Independent Normal proposal using fitted mean and covariance matrix
#'@export
get_proposal_independent_normal <- function(){
  f = function(thetas,normw,...){
    covariance = cov.wt(t(thetas), wt = normw, method = "ML")
    mean_t = covariance$center
    cov_t = covariance$cov + diag(rep(10^(-4)/nrow(thetas)), nrow(thetas)) # increased a bit the diagonal to prevent degeneracy effects)
    # define the sampler
    rproposal = function(Ntheta) {
      return (fast_rmvnorm_transpose(Ntheta, mean_t, cov_t))
    }
    # define the corresponding density function
    dproposal = function(thetas,log = TRUE) {
      if (log) {return (fast_dmvnorm_transpose(thetas, mean_t, cov_t))}
      else {return (exp(fast_dmvnorm_transpose(thetas, mean_t, cov_t)))}
    }
    return (list(r = rproposal, d = dproposal))
  }
  return(f)
}

#'@rdname get_proposal_mixture
#'@title get_proposal_mixture
#'@description Independent proposal from a fitted mixture of Normals with \code{nclust} components (default is 5).
#' If the fit is unsuccessful, return independent Normal proposal (see \code{get_independent_normal_proposal}).
#'@export
get_proposal_mixture <- function(nclust = 5, maxattempts = 5, verbose = FALSE){
  f <- function(thetas,normw,...){
    options(warn = -1)
    # resample
    ancestors <- systematic_resampling_n(normw, length(normw), runif(1))
    thetas_check <- thetas[,ancestors,drop=FALSE]
    # fit mixture
    fit <- mixmodCluster(data = data.frame(t(thetas_check)), nbCluster = nclust, dataType = "quantitative")
    # test that it worked
    is.error <- (length(fit@bestResult@parameters@proportions) == 0)
    attempt = 0
    while(attempt < maxattempts && is.error){
      attempt <- attempt + 1
      if(verbose){cat("fitting mixture... attempt", attempt, "\n")}
      fit <- mixmodCluster(data = data.frame(t(thetas_check)), nbCluster = nclust, dataType = "quantitative")
      # test that it worked
      is.error <- (length(fit@bestResult@parameters@proportions) == 0)
    }
    options(warn = 0)
    if (is.error){
      return(get_proposal_independent_normal()(thetas,normw))
    }

    # if it worked, ...
    rproposal = function(Ntheta) {
      proportions <- fit@bestResult@parameters@proportions
      means <- fit@bestResult@parameters@mean
      variances <- fit@bestResult@parameters@variance
      K <- nrow(means)
      X <- matrix(0, ncol = Ntheta, nrow = ncol(means))
      # sample allocations
      allocations <- systematic_resampling_n(proportions, Ntheta, runif(1))
      for (k in 1:K){
        which.k <- which(allocations == k)
        nk <- length(which.k)
        if (nk > 0){
          xk <- fast_rmvnorm_transpose(nk, means[k,], variances[[k]])
          X[,which.k] <- xk
        }
      }
      # random shuffling
      X <- X[,sample(x = 1:ncol(X), size = ncol(X), replace = FALSE),drop=FALSE]
      return(X)
    }
    dproposal = function(thetas,log = TRUE) {
      proportions <- fit@bestResult@parameters@proportions
      means <- fit@bestResult@parameters@mean
      variances <- fit@bestResult@parameters@variance
      d <- nrow(thetas)
      n <- ncol(thetas)
      K <- nrow(means)
      evals <- matrix(0, nrow = n, ncol = K)
      for (k in 1:K){
        evals[,k] <- fast_dmvnorm_transpose(thetas, means[k,], variances[[k]]) + log(proportions[k])
      }
      g <- function(row){
        m <- max(row)
        return(m + log(sum(exp(row - m))))
      }
      results <- apply(X = evals, MARGIN = 1, FUN = g)
      if (log) return(results)
      else return(exp(results))
    }
    return (list(r = rproposal, d = dproposal))
  }
  return(f)
}

