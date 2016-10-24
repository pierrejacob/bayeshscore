rm(list = ls())
library(HyvarinenSSM)
library(doMC)
registerDoMC(cores = 4)
set.seed(17)

model.pz5 <- get_model_pz5()
model.pz5$theta = c(0.7,0.5,0.25,0.3)
nobservations <- 100
datafile <- paste0("~/Dropbox/Harvard/Grant/NSF DMS Fall 2016/pzdata.T", nobservations, ".RData")
## generate data
# sim = simulateData(model.pz, theta = c(0.7,0.5,0.25,0.3), nobservations)
# save(sim, file = datafile)
## load data
load(datafile)
X = sim$X
Y = sim$Y


observations <- matrix(Y, nrow = model.pz5$dimY)# observations in a matrix of dimensions dimy x nobservations

algorithmic_parameters <- list(Ntheta = 2^10, Nx = 2^9, nmoves = 3,
                               resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
                               progress = TRUE)

savefile <- paste0("~/Dropbox/Harvard/Grant/NSF DMS Fall 2016/hscore.pz5.T", nobservations,
                  "Nt", algorithmic_parameters$Ntheta, "Nx", algorithmic_parameters$Nx,  ".RData")

nrep <- 4
results.pz5 = foreach(i=1:nrep,.packages='HyvarinenSSM',.verbose = TRUE) %dorng% {
  hscore_continuous(observations, model.pz5, algorithmic_parameters)
}
save(results.pz5, file = savefile)

names(results.pz5[[1]])
hscores <- foreach (irep = 1:nrep, .combine = cbind) %dopar% {
  results.pz5[[irep]]$hscore
}
evidences <- foreach (irep = 1:nrep, .combine = cbind) %dopar% {
  results.pz5[[irep]]$logevidence
}

ESSs <- foreach (irep = 1:nrep, .combine = cbind) %dopar% {
  results.pz5[[irep]]$ESS
}

matplot(hscores, type = "l")
matplot(evidences, type = "l")
matplot(ESSs, type = "l")
length(unique(results.pz5[[1]]$thetas))
