#'@rdname hscore_discrete
#'@title hscore_discrete
#'@description This function computes successive prequential Hyvarinen score for discrete observations by running the smc^2 algorithm. It also computes the successive log-evidence as a by-product.
#'@export
hscore_discrete <- function(observations, model, algorithmic_parameters){
  Nx <- algorithmic_parameters$Nx
  Ntheta <- algorithmic_parameters$Ntheta
  resampling <- algorithmic_parameters$resampling
  nobservations <- ncol(observations)
  if (algorithmic_parameters$progress) {
    progbar = txtProgressBar(min = 0,max = nobservations,style=3)
    count = 0
    time_start = proc.time()
  }
  Xpred = array(NA,dim = c(Nx, model$dimX, Ntheta)) #particles X w.r.t. predictive distributions (most recent)
  ESS = array(NA,dim = c(nobservations)) #ESS at successive times t
  Hscore = array(NA,dim = c(nobservations)) #prequential Hyvarinen score at successive times t
  logevidence = array(NA,dim = c(nobservations)) #log-evidence at successive times t
  # sample from prior
  thetas <- model$rprior(Ntheta)
  # initialize filter
  first_step <- filter_init(observations[,1], model, thetas, algorithmic_parameters)
  X = first_step$X  #Nx particles of dimension 1 for each theta (most recent)
  xnormW = first_step$xnormW  #matrix of normalized weights for X (most recent)
  Xpred = X
  XprednormW = matrix(1/Nx, nrow = Nx, ncol = Ntheta) #matrix of normalized weights for Xpred (most recent)
  z = first_step$z  #matrix of likelihood estimates
  thetalogw <- z #update log-weights for theta
  maxlogW <- max(thetalogw) #avoids overflow when exponentiating
  W <- exp(thetalogw - maxlogW) #computes actual unnormalized weights for theta
  logevidence[1] <- log(mean(W)) + maxlogW #udpate log-evidence
  thetanormw <- W / sum(W) #normalize weights for theta
  # compute H score
  Hscore[1] = SHd(1,model,observations[,1],thetas,rep(1/Ntheta,Ntheta),Xpred,XprednormW,Ntheta,Nx)
  # ...
  ESS[1] <- getESS(thetanormw)
  if ((ESS[1]/Ntheta) < 0.5){
    rejuvenation = rejuvenation_step(observations, 1, model, thetas, thetanormw, X, xnormW, z, algorithmic_parameters)
    thetas = rejuvenation$thetas
    thetalogw = rep(log(1/Ntheta),Ntheta)
    X = rejuvenation$X
    xnormW = rejuvenation$xnormW
    z = rejuvenation$z
  }
  if (algorithmic_parameters$progress) {
    count = count + 1
    setTxtProgressBar(progbar, count)
  }
  for (t in 2:nobservations){
    Xpred = filter_predict(t, model, thetas, X, algorithmic_parameters)
    XprednormW = xnormW
    # compute H score
    Hscore[t] = Hscore[t-1] + SHd(t,model,observations[,t],thetas,thetanormw,Xpred,XprednormW,Ntheta,Nx)
    #...
    next_step <- filter_next_step(observations[,t], t, model, thetas, X, xnormW, algorithmic_parameters)
    X = next_step$X  #Nx particles of dimension 1 for each theta (most recent)
    xnormW = next_step$xnormW  #matrix of normalized weights for X (most recent)
    z_incremental = next_step$z_incremental  #matrix of likelihood estimates
    z <- z + z_incremental
    # update log-evidence estimator
    maxz_incremental = max(z_incremental)
    logevidence[t] <- logevidence[t-1] + log(sum(exp(z_incremental - maxz_incremental)*thetanormw)) + maxz_incremental
    #
    thetalogw <- thetalogw + z_incremental #update log-weights for theta
    maxlogW <- max(thetalogw) #avoids overflow when exponentiating
    W <- exp(thetalogw - maxlogW) #computes actual unnormalized weights for theta
    thetanormw <- W / sum(W) #normalize weights for theta
    # ...
    ESS[t] <- getESS(thetanormw)
    if ((ESS[t]/Ntheta) < 0.5){
      rejuvenation = rejuvenation_step(observations, t, model, thetas, thetanormw, X, xnormW, z, algorithmic_parameters)
      thetas = rejuvenation$thetas
      thetalogw = rep(log(1/Ntheta),Ntheta)
      X = rejuvenation$X
      xnormW = rejuvenation$xnormW
      z = rejuvenation$z
    }
    if (algorithmic_parameters$progress) {
      count = count + 1
      setTxtProgressBar(progbar, count)
    }
  }
  if (algorithmic_parameters$progress) {
    close(progbar)
    time_end = proc.time()-time_start
    cat(paste("Hscore: T = ",toString(nobservations),", Ntheta = ",toString(Ntheta),
              ", Nx = ",toString(Nx),"\n",sep = ""))
    print(time_end)
  }
  return (list(hscore = Hscore, logevidence = logevidence, ESS = ESS, thetas = thetas, thetanormw = thetanormw))
}

