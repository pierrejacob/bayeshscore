rm(list = ls())
library(HyvarinenSSM)
library(doMC)
set.seed(17)

nobservations <- 50
model <- get_model_lineargaussian()
sim = simulateData(model, theta = c(0.8,1,1,1), nobservations)
X = sim$X
Y = sim$Y

# observations in a matrix of dimensions dimy x nobservations
observations <- matrix(Y, nrow = model$dimY)
algorithmc_parameters <- list(Nx = 1024, resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)))
theta <- model$theta

rhos <- seq(from = 0, to = 1, length.out = 100)
registerDoMC(cores = 10)
loglikelihoods <- foreach(rho = rhos, .combine = c) %dorng% {
  bpfresults <- bootstrap_particle_filter(observations, model, c(rho, theta[2:4]), algorithmc_parameters)
  bpfresults$log_p_y_hat
}

plot(x = rhos, y = loglikelihoods)


