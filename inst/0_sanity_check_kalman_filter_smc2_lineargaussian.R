rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(19)

# create data
nobservations <- 30
model <- get_model_simplerlineargaussian()
theta_star <- c(0.8,1,model$psi,model$sigmaV2)
sim = simulateData(model, theta = theta_star, nobservations)
X = sim$X
Y = sim$Y
observations <- matrix(Y, nrow = model$dimY)
# observations in a matrix of dimensions dimY by nobservations

# set algorithmic parameters
algorithmic_parameters <- list()
algorithmic_parameters$Ntheta = 1024
algorithmic_parameters$Nx = 128
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
model$dpredictive = NULL # this forces the use of SMC2
smc2_results <- hscore(observations, model, algorithmic_parameters)

# Plot the filtering means (and some particles X for one theta, e.g. the one with largest weight)
xmean = rep(NA,nobservations)
xmeanMAP = rep(NA,nobservations)
sim = data.frame()
indexMAP = which.max(smc2_results$normw_history[[nobservations+1]])
Ntheta = algorithmic_parameters$Ntheta
for (t in 1:nobservations){
  PF = smc2_results$PF_history[[t+1]]
  normw = smc2_results$normw_history[[t+1]]
  xmean[t] = sum(sapply(1:Ntheta,function(i)sum(PF[[i]]$X*PF[[i]]$xnormW)*normw[i]))
  xmeanMAP[t] = sum(PF[[indexMAP]]$X*PF[[indexMAP]]$xnormW)
  sim = rbind(sim, data.frame(X = c(PF[[indexMAP]]$X), w = c(PF[[indexMAP]]$xnormW), time = t))
}

# set parameters
thetaMAP = smc2_results$thetas_history[[nobservations+1]][,indexMAP]
phi = thetaMAP[1]
sigmaW2 = thetaMAP[2]
psi = thetaMAP[3]
sigmaV2 = thetaMAP[4]
initial_mean = 0
initial_var = (sigmaW2)/(1-phi^2)
# run Kalman filter
KF = KF_filtering(Y,phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var)


ggplot() +
  geom_point(data=sim,aes(time,X,alpha=w)) +
  # geom_line(aes(1:nobservations,X[,]),col='blue',linetype='dashed') +
  # geom_line(aes(1:nobservations,Y[,]),col='blue',size=1) +
  geom_line(aes(1:nobservations,sapply(1:nobservations,function(t)KF[[t]]$muX_t_t)),col='black',size=1) +
  geom_line(aes(1:nobservations,xmeanMAP),col='purple',size=2,linetype='dashed') +
  theme(legend.position="none") +
  ylab('particles X (dots), X (blue), Y (red)')
