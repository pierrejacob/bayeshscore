rm(list = ls())
library(HyvarinenSSM)
set.seed(17)

nobservations <- 100
model <- get_model_lineargaussian()
sim = simulateData(model, theta = c(0.8,1,1,1), nobservations)
X = sim$X
Y = sim$Y
observations <- matrix(Y, nrow = model$dimY)# observations in a matrix of dimensions dimy x nobservations


# Plot data
observations.df = data.frame(time = 1:nobservations, X = t(X),Y = t(Y))
g = ggplot(observations.df, aes(x = time)) +
  geom_point(aes(y=Y),size=2) +
  geom_line(aes(y=X),linetype=2) +
  xlab("\n Time")
plot(g)


algorithmic_parameters <- list(Ntheta = 2^10, Nx = 2^10,
                              resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
                              progress = TRUE)

# # define proposals for PMMH
# algorithmic_parameters$MHsd = 0.05
# algorithmic_parameters$rMHproposal = function(current_theta) {
#   return (fast_rmvnorm(1, current_theta, diag(algorithmic_parameters$MHsd^2,4)))
# }
# algorithmic_parameters$dMHproposal = function(new_theta,current_theta) {
#   #Note: this outputs the LOG-density
#   return (fast_dmvnorm(new_theta, current_theta, diag(algorithmic_parameters$MHsd^2,4)))
# }
# algorithmic_parameters$M = 1000 #number of initial PMMH iterations
# algorithmic_parameters$burnin = 500 #burn-in
# algorithmic_parameters$initialbatchsize = 1 #number of observations to include in initial PMMH

# initial_proposal = get_batch_initial_proposal(observations,model,algorithmic_parameters)
# logevidence_1_to_b = get_batch_logevidence(observations, model, initial_proposal, algorithmic_parameters)
# algorithmic_parameters$rinitial_theta = initial_proposal$r
# algorithmic_parameters$dinitial_theta = initial_proposal$logd

result = hscore_continuous(observations, model, algorithmic_parameters)

#sanity check posterior of theta (first component)
setmytheme()
qplot(x = result$thetas[,1], weight = result$thetanormw, geom = "histogram",binwidth = 0.05,xlim = c(0,1))

#sanity check
par(mfrow=c(3,1))
plot(result$ESS,type='l')
plot(result$logevidence,type='l')
plot(result$hscore,type='l')
par(mfrow=c(1,1))

print(result$logevidence)
print(result$hscore)
