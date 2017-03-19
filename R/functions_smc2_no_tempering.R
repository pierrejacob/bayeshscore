filter_first_step <- function(observation1, model, thetas, trees, algorithmic_parameters){
  # Extract algorithmic parameters
  Ntheta <- algorithmic_parameters$Ntheta
  Nx <- algorithmic_parameters$Nx
  # Initialize array of particles X and associated weights
  X = array(NA,dim = c(Nx, model$dimX, Ntheta)) #Nx particles (most recent) for each theta (size = Nx,dimX,Ntheta)
  xnormW = matrix(0, nrow = Nx, ncol = Ntheta) #matrix of corresponding normalized X-weights (size = Nx,Ntheta)
  log_z = rep(0, Ntheta) #matrix of log-likelihood estimates (size = Ntheta)
  for (itheta in 1:Ntheta){
    X[,,itheta] <- model$rinitial(thetas[itheta,],Nx) #initial step 1
    logW <- model$dobs(observation1, X[,,itheta], 1, thetas[itheta,],log = TRUE)
    maxlogW <- max(logW)
    W <- exp(logW - maxlogW)
    log_z[itheta] <- log(mean(W)) + maxlogW #udpate likelihood estimate
    xnormW[,itheta] <- W / sum(W)
    # Store paths in tree
    trees[[itheta]]$init(t(X[,,itheta]))
  }
  return(list(X = X, xnormW = xnormW, log_z = log_z, trees = trees))
}

filter_next_step <- function(observationt, t, model, thetas, X, xnormW, trees, algorithmic_parameters){
  Ntheta <- algorithmic_parameters$Ntheta
  resampling <- algorithmic_parameters$resampling
  Xnew = X
  xnormWnew = xnormW
  log_z_incremental = rep(0, Ntheta) #matrix of likelihood estimates
  for (itheta in 1:Ntheta){
    ancestors <- resampling(xnormW[,itheta]) #sample the ancestors' indexes
    X_current <- X[ancestors,,itheta]
    if (is.null(dim(X_current))) X_current <- matrix(X_current, ncol = model$dimX)
    X_current <- model$rtransition(X_current, t, thetas[itheta,])
    logW <- model$dobs(observationt, X_current, t, thetas[itheta,],log = TRUE)
    maxlogW <- max(logW)
    W <- exp(logW - maxlogW)
    log_z_incremental[itheta] <- log(mean(W)) + maxlogW #udpate likelihood estimate
    xnormWnew[,itheta] <- W / sum(W)
    Xnew[,,itheta] <- X_current
    trees[[itheta]]$update(t(Xnew[,,itheta]), ancestors-1)
  }
  return(list(X = Xnew, xnormW = xnormWnew, log_z_incremental = log_z_incremental, trees = trees))
}

rejuvenation_step <- function(observations, t, model, thetas, thetanormw, X, xnormW, log_z, trees, algorithmic_parameters){
  # Extract algorithmic parameters
  Ntheta = algorithmic_parameters$Ntheta
  Nx = algorithmic_parameters$Nx
  nmoves = algorithmic_parameters$nmoves
  if (is.null(nmoves)) nmoves = 1
  resampling = algorithmic_parameters$resampling
  # Compute parameters for proposal move step
  covariance = cov.wt(thetas,wt = thetanormw,method = "ML")
  mean_t = covariance$center
  cov_t = matrix(covariance$cov,nrow = model$dimtheta) + diag(rep(10^(-4)/model$dimtheta),model$dimtheta)
  #(increased a little bit the diagonal to prevent degeneracy effects)
  resampled_index = resampling(thetanormw)
  thetas = thetas[resampled_index,,drop=FALSE]
  X = X[,,resampled_index,drop=FALSE]
  xnormW = xnormW[,resampled_index]
  log_z = log_z[resampled_index]
  trees = trees[resampled_index]
  #
  for (imove in 1:nmoves){
    theta_new_all = fast_rmvnorm(Ntheta,mean_t,cov_t)
    proposal_density_new_all <- fast_dmvnorm(theta_new_all, mean_t, cov_t)
    proposal_density_current <- fast_dmvnorm(thetas, mean_t, cov_t)
    accepts <- 0
    for (i in 1:Ntheta) {
      theta_old <- thetas[i,]
      theta_new <- theta_new_all[i,]
      log_z_old = log_z[i]
      logprior_theta_old <- model$dprior(theta_old, log = TRUE) # wasteful; this has been computed before ...
      logprior_theta_new <- model$dprior(theta_new, log = TRUE)
      if (logprior_theta_new == -Inf){
        next
      }
      else {
        #the CPF function performs a regular PF when no conditioning path is provided
        PF <- conditional_particle_filter(matrix(observations[,1:t],ncol = t), model, theta_new, Nx)
        log_z_new <- PF$log_p_y_hat
        lognum <- logprior_theta_new + log_z_new + proposal_density_current[i]
        logdenom <- logprior_theta_old + log_z_old + proposal_density_new_all[i]
        logacceptance <- lognum - logdenom
        logu <- log(runif(1))
        if (logu <= logacceptance){
          accepts <- accepts + 1
          thetas[i,] <- theta_new
          X[,,i] <- PF$X
          xnormW[,i] <- PF$xnormW
          log_z[i] <- log_z_new
          trees[[i]] = PF$tree
        }
        else {
          # do nothing
        }
      }
    }
  }
  return(list(thetas = thetas, X = X, xnormW = xnormW, log_z = log_z, trees = trees, accept_rate = accepts/Ntheta))
}

filter_predict <- function(t, model, thetas, X, algorithmic_parameters){
  Ntheta <- algorithmic_parameters$Ntheta
  Nx <- algorithmic_parameters$Nx
  Xpred = X
  for (itheta in 1:Ntheta){
    if (is.null(dim(X[,,itheta]))){
      Xpred[,,itheta] = model$rtransition(matrix(X[,,itheta],nrow = Nx), t, thetas[itheta,])
    }
    else{
      Xpred[,,itheta] = model$rtransition(X[,,itheta], t, thetas[itheta,])
    }
  }
  return (Xpred)
}
