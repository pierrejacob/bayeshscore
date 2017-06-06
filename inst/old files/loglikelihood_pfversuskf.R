rm(list = ls())
library(HyvarinenSSM)
library(doMC)
set.seed(17)
set_global_path()

# # load data
# load(file = paste0(rdatapath, "lineargaussian.observations.RData"))
# #
# nobservations <- 50
# X <- X[,1:nobservations,drop=F]
# Y <- Y[,1:nobservations,drop=F]
# # load model
# model <- get_model_lineargaussian()
# # observations in a matrix of dimensions dimy x nobservations
# observations <- matrix(Y, nrow = model$dimY)


nobservations <- 50
model <- get_model_simplerlineargaussian()
theta_star <- model$theta
sim = simulateData(model, theta = theta_star, nobservations)
X = sim$X
Y = sim$Y
# observations in a matrix of dimensions dimy x nobservations
observations <- matrix(Y, nrow = model$dimY)
algorithmc_parameters <- list(Nx = 2^14, resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)))
theta <- theta_star

# calculate log-likelihood for various parameters (changing only the first component)
library(dlm)
phi <- theta[1]
psi = model$psi
sigmaV2 = model$sigmaV2
sigmaW2 <- theta[2]

mod <- dlm(list(m0 = 0, C0 = sigmaW2/(1-phi^2), FF = psi, V = sigmaV2, GG = phi, W = sigmaW2))
-dlmLL(observations[1,], mod)
library(astsa)
astsa_results <- Kfilter0(nobservations, observations[1,], psi, 0, sigmaW2/(1-phi^2), phi, sqrt(sigmaW2), sqrt(sigmaV2))
-astsa_results$like

kalman_module <- Module( "kalman_mod", PACKAGE = "HyvarinenSSM")
Kalman <- new(kalman_module$Kalman)
Kalman$set_parameters(list(rho = theta[1], sigma = sqrt(theta[4]), eta = theta[2], tau = sqrt(theta[3])))
Kalman$set_observations(matrix(observations, ncol = 1))
Kalman$filtering()
Kalman$getLL()
Kalman$get_filtering_mean(10)
Kalman$get_filtering_mean(4)
astsa_results$xf[,,10]
astsa_results$xf[,,4]

kalman_ll <- function(theta){
  phi <- theta[1]
  psi = model$psi
  sigmaV2 = model$sigmaV2
  sigmaW2 <- theta[2]
  Kalman <- new(kalman_module$Kalman)
  Kalman$set_parameters(list(rho = phi, sigma = sqrt(sigmaW2), eta = psi, tau = sqrt(sigmaV2)))
  Kalman$set_observations(matrix(observations, ncol = 1))
  Kalman$first_step()
  for (i in 1:nobservations){
    Kalman$filtering_step(i-1)
  }
  return(Kalman$get_incremental_ll())
}

kalman_ll(theta_star)


nparam <- 10
theta1grid <- seq(from = 0.3, to = 0.9, length.out = nparam)
log_p_y_hat <- c()
log_p_y <- c()
log_p_y2 <- c()
for (i in 1:nparam){
  theta[1] <- theta1grid[i]
  #
  log_p_y <- c(log_p_y, sum(kalman_ll(theta)))
  #
  pf_results <- bootstrap_particle_filter(observations, model, theta, algorithmc_parameters)
  log_p_y_hat <- c(log_p_y_hat, pf_results$log_p_y_hat)
  phi <- theta[1]
  psi = model$psi
  sigmaV2 = model$sigmaV2
  sigmaW2 <- theta[2]
  astsa_results <- Kfilter0(nobservations, observations[1,], psi, 0, sigmaW2/(1-phi^2), phi, sqrt(sigmaW2), sqrt(sigmaV2))
  log_p_y2 <- c(log_p_y2, -astsa_results$like)

}
plot(theta1grid, log_p_y_hat)
lines(theta1grid, log_p_y)


log_p_y - log_p_y_hat
log_p_y - log_p_y2
log_p_y_hat - log_p_y2
