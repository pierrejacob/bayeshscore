rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(19)

# create data
nobservations <- 50
model <- get_model_simplerlineargaussian()
theta_star <- c(0.8,1)
sim = simulateData(model, theta = theta_star, nobservations)
X = sim$X
Y = sim$Y
observations <- matrix(Y, nrow = model$dimY)
# observations in a matrix of dimensions dimY by nobservations

# set algorithmic parameters
algorithmic_parameters <- list()
algorithmic_parameters$Ntheta = 1024
algorithmic_parameters$Nx = 128
algorithmic_parameters$observation_type = 'continuous'
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$store_X = TRUE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$min_acceptance_rate = 0.45
algorithmic_parameters$nmoves = 2
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R

### Run SMC_2
module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree
model$likelihood = NULL # this forces the use of SMC2
smc2_results <- hscore(observations, model, model)


# index = algorithmic_parameters$Ntheta
xmean = rep(NA,nobservations)
sim = data.frame()
indexMAP = which.max(smc2_results_temp$normw_history[[nobservations+1]])
for (t in 1:nobservations){
  PF = smc2_results_temp$PF_history[[t+1]]
  xmean[t] = 0
  for (i in 1:algorithmic_parameters$Ntheta){
    xmean[t] = xmean[t] + sum(PF[[i]]$X*PF[[i]]$xnormW)*smc2_results_temp$normw_history[[t+1]][i]
  }
  sim = rbind(sim, data.frame(X = c(PF[[indexMAP]]$X),
                              w = c(PF[[indexMAP]]$xnormW),
                              time = t))
}
ggplot() +
  geom_point(data=sim,aes(time,X,alpha=w)) +
  geom_line(aes(1:nobservations,c(X)),col='blue',size=2) +
  geom_line(aes(1:nobservations,c(Y)),col='red',size=2) +
  geom_line(aes(1:nobservations,xmean),col='black',size=2) +
  theme(legend.position="none") +
  ylab('particles X (dots), X (blue), Y (red)')
