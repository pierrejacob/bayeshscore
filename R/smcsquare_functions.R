filter_init <- function(observation1, model, thetas, algorithmic_parameters){
  Ntheta <- algorithmic_parameters$Ntheta
  Nx <- algorithmic_parameters$Nx
  X = array(NA,dim = c(Nx, model$dimX, Ntheta)) #Nx particles of dimension dimX for each theta (most recent)
  z = rep(0, Ntheta) #matrix of likelihood estimates
  xnormW = matrix(NA, nrow = Nx, ncol = Ntheta) #matrix of normalized weights for X (most recent)
  for (itheta in 1:Ntheta){
    X[,,itheta] <- model$rinitial(thetas[itheta,],Nx) #initial step 1
    logW <- model$dobs(observation1, X[,,itheta], 1, thetas[itheta,],log = TRUE)
    maxlogW <- max(logW)
    W <- exp(logW - maxlogW)
    z[itheta] <- log(mean(W)) + maxlogW #udpate likelihood estimate
    xnormW[,itheta] <- W / sum(W)
  }
  return(list(X = X, xnormW = xnormW, z = z))
}

filter_next_step <- function(observationt, t, model, thetas, X, xnormW, algorithmic_parameters){
  Ntheta <- algorithmic_parameters$Ntheta
  resampling <- algorithmic_parameters$resampling
  Xnew = X
  xnormWnew = xnormW
  z_incremental = rep(0, Ntheta) #matrix of likelihood estimates
  for (itheta in 1:Ntheta){
    ancestors <- resampling(xnormW[,itheta]) #sample the ancestors' indexes
    X_current <- X[ancestors,,itheta]
    if (is.null(dim(X_current))) X_current <- matrix(X_current, ncol = model$dimX)
    X_current <- model$rtransition(X_current, t, thetas[itheta,])
    logW <- model$dobs(observationt, X_current, t, thetas[itheta,],log = TRUE)
    maxlogW <- max(logW)
    W <- exp(logW - maxlogW)
    z_incremental[itheta] <- log(mean(W)) + maxlogW #udpate likelihood estimate
    xnormWnew[,itheta] <- W / sum(W)
    Xnew[,,itheta] <- X_current
  }
  return(list(X = Xnew, xnormW = xnormWnew, z_incremental = z_incremental))
}

rejuvenation_step <- function(observations, t, model, thetas, thetanormw, X, xnormW, z, algorithmic_parameters){
  Ntheta <- algorithmic_parameters$Ntheta
  Nx <- algorithmic_parameters$Nx
  resampling <- algorithmic_parameters$resampling
  thetasnew <- matrix(NA, nrow = Ntheta, ncol = model$dimtheta)
  Xnew <- array(NA,dim = dim(X))
  xnormWnew <-  matrix(NA, nrow = Nx, ncol = Ntheta)
  znew = rep(0, Ntheta)
  #compute parameters for proposal move step
  covariance = cov.wt(thetas,wt = thetanormw)
  mean_t <- covariance$center
  cov_t <- matrix(covariance$cov,nrow = model$dimtheta)
  # increased a little bit the diagonal to prevent degeneracy effects
  # cov_t <- matrix(covariance$cov,nrow = model$dimtheta) + diag(rep(10^(-4)/model$dimtheta),model$dimtheta)
  resampled_index = resampling(thetanormw)
  theta_new_all = fast_rmvnorm(Ntheta,mean_t,cov_t)
  for (i in 1:Ntheta) {
    theta_old <- thetas[resampled_index[i],]
    theta_new <- theta_new_all[i,]
    z_old = z[resampled_index[i]]
    logprior_theta_old <- model$dprior(theta_old, log = TRUE)
    logprior_theta_new <- model$dprior(theta_new, log = TRUE)
    if (logprior_theta_new == -Inf){
      thetasnew[i,] <- theta_old
      Xnew[,,i] <- X[,,resampled_index[i]]
      xnormWnew[,i] <- xnormW[,resampled_index[i]]
      znew[i] <- z_old
      next
    }
    PF <- bootstrap_particle_filter(matrix(observations[,1:t],ncol = t), model, theta_new, algorithmic_parameters)
    z_new <- PF$log_p_y_hat
    lognum <- logprior_theta_new + z_new + fast_dmvnorm(matrix(theta_old,ncol = model$dimtheta),mean_t,cov_t)
    logdenom <- logprior_theta_old + z_old + fast_dmvnorm(matrix(theta_new,ncol = model$dimtheta),mean_t,cov_t)
    logacceptance <- lognum - logdenom
    logu <- log(runif(1))
    if (logu <= logacceptance){
      thetasnew[i,] <- theta_new
      Xnew[,,i] <- PF$X
      xnormWnew[,i] <- PF$normW
      znew[i] <- z_new
    }
    else {
      thetasnew[i,] <- theta_old
      Xnew[,,i] <- X[,,resampled_index[i]]
      xnormWnew[,i] <- xnormW[,resampled_index[i]]
      znew[i] <- z_old
    }
  }
  return(list(thetas = thetasnew, X = Xnew, xnormW = xnormWnew, z = znew))
}
