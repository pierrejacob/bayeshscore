rm(list = ls())
library(HyvarinenSSM)
library(doMC)
set.seed(17)

nobservations <- 25
model <- get_model_lineargaussian()
sim = simulateData(model, theta = c(0.8,1,1,1), nobservations)
X = sim$X
Y = sim$Y

# observations in a matrix of dimensions dimy x nobservations
observations <- matrix(Y, nrow = model$dimY)
algorithmc_parameters <- list(Nx = 2^10, resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)))
theta <- model$theta

module_tree <- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <- module_tree$Tree

conditional_particle_filter <- function(observations, model, theta, algorithmic_parameters, path = NULL){
  Tree <- new(TreeClass, algorithmc_parameters$Nx, 10*algorithmc_parameters$Nx, model$dimX)
  # we expect a path as a matrix of size dimX x nobservations
  cpf <- TRUE
  if (is.null(path)){
    cpf <- FALSE
  }
  Nx <- algorithmic_parameters$Nx
  # resampling <- algorithmic_parameters$resampling
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
  if (cpf){
    X[Nx,] <- path[,1]
  }
  Tree$init(t(X))
  #iterate for n = 2, ... T
  if (nobservations > 1){
    for (t in 2:nobservations) {
      ancestors <- sample(x = 1:Nx, size = Nx, replace = TRUE, prob = normW) #sample the ancestors' indexes
      if (cpf){
        ancestors[Nx] <- Nx
      }
      X <- X[ancestors,]
      if (is.null(dim(X))) X <- matrix(X, ncol = model$dimX)
      X <- model$rtransition(X, t, theta)
      if (cpf){
        X[Nx,] <- path[,t]
      }
      logW <- model$dobs(observations[,t], X, t, theta)
      maxlogW <- max(logW)
      W <- exp(logW - maxlogW)
      log_p_y_hat <- log_p_y_hat + log(mean(W)) + maxlogW #udpate likelihood estimate
      normW <- W / sum(W)
      Tree$update(t(X), ancestors-1)
    }
  }
  new_path <- Tree$get_path(sample(x = 0:(Nx-1), size = 1, replace = TRUE, prob = normW))
  return(list(log_p_y_hat = log_p_y_hat, X = X, normW = normW, path = new_path))
}

cpf_results <- conditional_particle_filter(observations, model, theta, algorithmc_parameters)
path <- cpf_results$path

niterations <- 1000
paths <- matrix(nrow = niterations, ncol = nobservations)
for (iteration in 1:niterations){
  cpf_results <- conditional_particle_filter(observations, model, theta, algorithmc_parameters, path = path)
  path <- cpf_results$path
  paths[iteration,] <- path[1,]
}

# matplot(paths[,1:10], type = "l")
cpf_smoothingmeans <- apply(paths, 2, mean)

library(dlm)
phi <- theta[1]
psi = theta[2]
sigmaV2 = theta[3]
sigmaW2 <- theta[4]

mod <- dlm(list(m0 = 0, C0 = sigmaW2/(1-phi^2), FF = psi, V = sigmaV2, GG = phi, W = sigmaW2))
dlmfilterresults <- dlmFilter(observations[1,], mod)
dlmsmootherresults <- dlmSmooth(observations[1,], mod)

plot(cpf_smoothingmeans, col = "blue")
lines(dlmsmootherresults$s[2:(nobservations+1)], col = "red")

