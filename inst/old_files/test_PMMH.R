rm(list = ls())
library(HyvarinenSSM)
library(gridExtra)
set.seed(17)

# Define model and data
observations <- data_kangaroo[c("y1","y2"),]
rangeprior = 10
model1 <- get_model_kangarooLogistic(rangeprior)
model2 <- get_model_kangarooExponential(rangeprior)
model3 <- get_model_kangarooRandomwalk(rangeprior)

# Define initial proposal for theta (to avoid sampling from vague prior)

# Define algorithmic parameters for each model
algorithmic_parameters = list(Ntheta = 2^10, Nx = 2^10,
              resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
              progress = TRUE)
# define proposals for PMMH
algorithmic_parameters$rMHproposal = function(current_theta) {
  return (fast_rmvnorm(1, current_theta, diag(c(0.05^2,0.005^2,0.05^2,0.0005^2))))
}
algorithmic_parameters$dMHproposal = function(new_theta,current_theta) {
  #Note: this outputs the LOG-density
  return (fast_dmvnorm(new_theta, current_theta, diag(c(0.1^2,0.01^2,0.1^2,0.001^2))))
}
algorithmic_parameters$M = 1000 #number of initial PMMH iterations
algorithmic_parameters$burnin = 500 #burn-in
algorithmic_parameters$initialbatchsize = 41 #number of observations to include in initial PMMH

PMMHsample = PMMH(observations,model1,algorithmic_parameters)

