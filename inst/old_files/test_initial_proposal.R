rm(list = ls())
library(HyvarinenSSM)
set.seed(17)

nobservations <- 50
model <- get_model_lineargaussian()
sim = simulateData(model, theta = c(0.8,1,1,1), nobservations)
X = sim$X
Y = sim$Y

# observations in a matrix of dimensions dimy x nobservations
observations <- matrix(Y, nrow = model$dimY)
algorithmic_parameters <- list(Nx = 1024, resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)))

# define proposals for PMMH
algorithmic_parameters$MHsd = 0.05
algorithmic_parameters$rMHproposal = function(current_theta) {
  return (fast_rmvnorm(1, current_theta, diag(algorithmic_parameters$MHsd^2,4)))
}
algorithmic_parameters$dMHproposal = function(new_theta,current_theta) {
  #Note: this outputs the LOG-density
  return (fast_dmvnorm(new_theta, current_theta, diag(algorithmic_parameters$MHsd^2,4)))
}
# set algorithmic parameters
algorithmic_parameters$M = 10000 #number of initial PMMH iterations
algorithmic_parameters$burnin = 5000 #burn-in
algorithmic_parameters$initialbatchsize = 1 #number of observations to include in initial PMMH
algorithmic_parameters$Ntheta = 2^10 #number of observations to include in initial PMMH

# sanity check PMMH
observations_1_to_b = matrix(observations[,1:algorithmic_parameters$initialbatchsize],nrow = model$dimY)
thetas = PMMH(observations_1_to_b,model,algorithmic_parameters)

# sanity check proposal and log-evidence
initial_proposal = get_batch_initial_proposal(observations,model,algorithmic_parameters)
logevidence_1_to_b = get_batch_logevidence(observations, model, initial_proposal, algorithmic_parameters)


