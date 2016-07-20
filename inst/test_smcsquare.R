rm(list = ls())
library(HyvarinenSSM)
set.seed(17)

nobservations <- 20
model <- get_model_lineargaussian()
sim = simulateData(model, theta = c(0.8,1,1,1), nobservations)
X = sim$X
Y = sim$Y

# observations in a matrix of dimensions dimy x nobservations
observations <- matrix(Y, nrow = model$dimY)
algorithmc_parameters <- list(Ntheta = 2^10, Nx = 1024,
                              resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)))

filter_init <- function(observation1, model, thetas, algorithmc_parameters){
  X = array(NA,dim = c(Nx, model$dimX, Ntheta)) #Nx particles of dimension 1 for each theta (most recent)
  z = rep(0, Ntheta) #matrix of likelihood estimates
  xnormW = matrix(0, nrow = Nx, ncol = Ntheta) #matrix of normalized weights for X (most recent)
  for (itheta in 1:algorithmc_parameters$Ntheta){
    X[,,itheta] <- model$rinitial(thetas[itheta,],Nx) #initial step 1
    logW <- model$dobs(observation1, X[,,itheta], 1, thetas[itheta,])
    maxlogW <- max(logW)
    W <- exp(logW - maxlogW)
    z[itheta] <- log(mean(W)) + maxlogW #udpate likelihood estimate
    xnormW[,itheta] <- W / sum(W)
  }
  return(list(X = X, xnormW = xnormW, z = z))
}

filter_next_step <- function(observationt, t, model, thetas, X, xnormW, algorithmc_parameters){
  Xnew = X
  xnormWnew = xnormW
  z_incremental = rep(0, Ntheta) #matrix of likelihood estimates
  # xnormW = matrix(nrow = Nx, ncol = Ntheta) #matrix of normalized weights for X (most recent)
  for (itheta in 1:algorithmc_parameters$Ntheta){
    ancestors <- resampling(xnormW[,itheta]) #sample the ancestors' indexes
    X_current <- X[ancestors,,itheta]
    if (is.null(dim(X_current))) X_current <- matrix(X_current, ncol = model$dimX)
    X_current <- model$rtransition(X_current, t, thetas[itheta,])
    logW <- model$dobs(observationt, X_current, t, thetas[itheta,])
    maxlogW <- max(logW)
    W <- exp(logW - maxlogW)
    z_incremental[itheta] <- log(mean(W)) + maxlogW #udpate likelihood estimate
    xnormWnew[,itheta] <- W / sum(W)
    Xnew[,,itheta] <- X_current
  }
  return(list(X = Xnew, xnormW = xnormWnew, z_incremental = z_incremental))
}

rejuvenation_step <- function(observations, t, model, thetas, X, xnormW, z, algorithmc_parameters){
  ....
  return(list(thetas = thetasnew, X = Xnew, xnormW = xnormWnew, z = newz))
}

## smcsquare <- function(...) ...

Nx <- algorithmc_parameters$Nx
Ntheta <- algorithmc_parameters$Ntheta
resampling <- algorithmc_parameters$resampling
theta <- model$theta

nobservations <- ncol(observations)

ESS = rep(0, nobservations) #ESS at successive times t
# sample from prior
thetas <- model$rprior(Ntheta)
# initialize filter
first_step <- filter_init(observations[,1], model, thetas, algorithmc_parameters)
X = first_step$X  #Nx particles of dimension 1 for each theta (most recent)
xnormW = first_step$xnormW  #matrix of normalized weights for X (most recent)
z = first_step$z  #matrix of likelihood estimates
#
thetalogw <- z
maxlogW <- max(thetalogw)
W <- exp(thetalogw - maxlogW)
thetanormw <- W / sum(W)
# if we want evidence estimator
evidence <- log(mean(W)) + maxlogW #udpate likelihood estimate
# compute H score here too
# ...
for (t in 2:nobservations){
  next_step <- filter_next_step(observations[,t], t, model, thetas, X, xnormW, algorithmc_parameters)
  X = next_step$X  #Nx particles of dimension 1 for each theta (most recent)
  xnormW = next_step$xnormW  #matrix of normalized weights for X (most recent)
  z_incremental = next_step$z_incremental  #matrix of likelihood estimates
  z <- z + z_incremental
  #
  thetalogw <- thetalogw + z_incremental
  maxlogW <- max(thetalogw)
  W <- exp(thetalogw - maxlogW)
  thetanormw <- W / sum(W)
  #
  # update evidence estimator
  # ...evidence <- evidence + ...
  # compute H score here too
  # ...
}

thetas
getESS(thetanormw)

setmytheme()
qplot(x = thetas[,1], weight = thetanormw, geom = "histogram")

