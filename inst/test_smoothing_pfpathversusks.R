rm(list = ls())
library(HyvarinenSSM)
library(doMC)
set.seed(17)

nobservations <- 10
model <- get_model_lineargaussian()
sim = simulateData(model, theta = c(0.8,1,1,1), nobservations)
X = sim$X
Y = sim$Y

# observations in a matrix of dimensions dimy x nobservations
observations <- matrix(Y, nrow = model$dimY)
algorithmc_parameters <- list(Nx = 2^14, resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)))
theta <- model$theta

module_tree <- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <- module_tree$Tree
Tree <- new(TreeClass, algorithmc_parameters$Nx, 10*algorithmc_parameters$Nx, model$dimX)

bootstrap_particle_filter <- function(observations, model, theta, algorithmic_parameters){
  Nx <- algorithmic_parameters$Nx
  resampling <- algorithmic_parameters$resampling
  nobservations <- ncol(observations)
  log_p_y_hat <- 0 #initialize estimate of p_theta(y_{1:nobservations}))
  X = matrix(NA,nrow = Nx, ncol = model$dimX) # matrix of Nx particles row-wise
  normW = rep(1/Nx, Nx) #vector of normalized weights
  X <- model$rinitial(theta,Nx) #initial step 1
  logW <- model$dobs(observations[,1],X,1,theta)
  maxlogW <- max(logW)
  W <- exp(logW - maxlogW)
  log_p_y_hat <- log_p_y_hat + log(mean(W)) + maxlogW #udpate likelihood estimate
  normW <- W / sum(W)
  #
  filtering_means <- rep(0, nobservations)
  filtering_means[1] <- sum(normW * X[,1])
  Tree$init(t(X))
  #iterate for n = 2, ... T
  if (nobservations > 1){
    for (t in 2:nobservations) {
      ancestors <- resampling(normW) #sample the ancestors' indexes
      X <- X[ancestors,]
      if (is.null(dim(X))) X <- matrix(X, ncol = model$dimX)
      X <- model$rtransition(X, t, theta)
      logW <- model$dobs(observations[,t], X, t, theta)
      maxlogW <- max(logW)
      W <- exp(logW - maxlogW)
      log_p_y_hat <- log_p_y_hat + log(mean(W)) + maxlogW #udpate likelihood estimate
      normW <- W / sum(W)
      filtering_means[t] <- sum(normW * X[,1])
      Tree$update(t(X), ancestors-1)
    }
  }
  return(list(log_p_y_hat = log_p_y_hat, X = X, normW = normW, filtering_means = filtering_means, Tree = Tree))
}

pf_results <- bootstrap_particle_filter(observations, model, theta, algorithmc_parameters)
plot(pf_results$filtering_means, type = "l")

library(dlm)
phi <- theta[1]
psi = theta[2]
sigmaV2 = theta[3]
sigmaW2 <- theta[4]

mod <- dlm(list(m0 = 0, C0 = sigmaW2/(1-phi^2), FF = psi, V = sigmaV2, GG = phi, W = sigmaW2))
dlmfilterresults <- dlmFilter(observations[1,], mod)
lines(dlmfilterresults$m[2:(nobservations+1)], col = "red")
dlmsmootherresults <- dlmSmooth(observations[1,], mod)

smoothing_means <- rep(0, nobservations)
for (iparticle in 1:algorithmc_parameters$Nx){
  trajectory <- pf_results$Tree$get_path(iparticle-1)
  smoothing_means <- smoothing_means + trajectory * pf_results$normW[iparticle]
}

lines(smoothing_means[1,], col = "blue")

plot(smoothing_means[1,], col = "blue")
lines(dlmsmootherresults$s[2:(nobservations+1)], col = "red")

