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
algorithmic_parameters <- list(Ntheta = 2^10, Nx = 2^10,
                              resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
                              progress = TRUE)
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


