rm(list = ls())
library(HyvarinenSSM)
library(doMC)
registerDoMC(cores = 4)
set.seed(17)

model.pz6 <- get_model_pz6()
model.pz6$theta = c(0.7,0.5,0.25,0.3)
nobservations <- 100
datafile <- paste0("~/Dropbox/Harvard/Grant/NSF DMS Fall 2016/pzdata.T", nobservations, ".RData")
## generate data
# sim = simulateData(model.pz, theta = c(0.7,0.5,0.25,0.3), nobservations)
# save(sim, file = datafile)
## load data
load(datafile)
X = sim$X
Y = sim$Y


observations <- matrix(Y, nrow = model.pz6$dimY)# observations in a matrix of dimensions dimy x nobservations

algorithmic_parameters <- list(Ntheta = 2^10, Nx = 2^9, nmoves = 3,
                               resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
                               progress = TRUE)

savefile <- paste0("~/Dropbox/Harvard/Grant/NSF DMS Fall 2016/hscore.pz6.T", nobservations,
                  "Nt", algorithmic_parameters$Ntheta, "Nx", algorithmic_parameters$Nx,  ".RData")

nrep <- 10
results.pz6 = foreach(i=1:nrep,.packages='HyvarinenSSM',.verbose = TRUE) %dorng% {
  hscore_continuous(observations, model.pz6, algorithmic_parameters)
}
save(results.pz6, file = savefile)

names(results.pz6[[1]])
hscores <- foreach (irep = 1:nrep, .combine = cbind) %dopar% {
  results.pz6[[irep]]$hscore
}
evidences <- foreach (irep = 1:nrep, .combine = cbind) %dopar% {
  results.pz6[[irep]]$logevidence
}

ESSs <- foreach (irep = 1:nrep, .combine = cbind) %dopar% {
  results.pz6[[irep]]$ESS
}

matplot(hscores, type = "l")
matplot(evidences, type = "l")
matplot(ESSs, type = "l")
length(unique(results.pz6[[1]]$thetas))
