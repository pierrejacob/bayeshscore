rm(list = ls())
library(HyvarinenSSM)
library(doMC)
set.seed(17)
set_global_path()

nobservations <- 1000
model <- get_model_lineargaussian()
theta_star <- model$theta
sim = simulateData(model, theta = theta_star, nobservations)
X = sim$X
Y = sim$Y

save(nobservations, theta_star, X, Y, file = paste0(rdatapath, "lineargaussian.observations.RData"))


