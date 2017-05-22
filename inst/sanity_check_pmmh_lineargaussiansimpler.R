##################################################################################################
# This checks that the outputs using SMC and SMC2 match the exact results
# in a linear gaussian case with 1 parameter (with discrete prior so that everything can be
# computed exactly).
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(19)

# create data
nobservations = 30
model = get_model_lineargaussiansimpler()
theta_star = c(0.8,1,model$psi,model$sigmaV2)
sim = simulateData(model, theta = theta_star, nobservations)
X = sim$X
Y = sim$Y
observations = matrix(Y, nrow = model$dimY)
# observations in a matrix of dimensions dimY by nobservations
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^14
algorithmic_parameters$Nx = 2^7
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$store_X = FALSE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$min_acceptance_rate = 0.0
algorithmic_parameters$nmoves = 1
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R

#--------------------------------------------------------------------------------------------
### Run SMC
smc_results = hscore(observations, model, algorithmic_parameters)
########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################
thetas_smc = smc_results$thetas_history[[nobservations+1]]
normw_smc = smc_results$normw_history[[nobservations+1]]
#--------------------------------------------------------------------------------------------
### Run PMMH
# define proposals for PMMH
PMMHcov = (cov.wt(t(thetas_smc[1:4,]),wt=smc_results$normw_history[[nobservations+1]])$cov)/5
PMMH_parameters = list()
PMMH_parameters$resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1))
PMMH_parameters$rMHproposal = function(current_theta) {
  return (fast_rmvnorm_transpose(1, current_theta, PMMHcov))
}
PMMH_parameters$dMHproposal = function(new_theta,current_theta) {
  #Note: this outputs the LOG-density
  return (fast_dmvnorm_transpose(new_theta, current_theta, PMMHcov))
}
PMMH_parameters$Nx = 100
PMMH_parameters$M = 3000 #number of initial PMMH iterations
PMMH_parameters$burnin = 1000 #burn-in
PMMH_parameters$progress = TRUE
thetasPMMH = PMMH(observations, model, PMMH_parameters)
#########################################################################################
#########################################################################################
#### Sanity check: sample posterior via naive MH
# Computes the posterior density (target)
psi = model$psi
sigmaV2 = model$sigmaV2
kalman_module <<- Module( "kalman_mod", PACKAGE = "HyvarinenSSM")
dpost = function(theta, log = TRUE){
  Kalman <<- new(kalman_module$Kalman)
  Kalman$set_parameters(list(rho = theta[1], sigma = sqrt(theta[2]), eta = psi, tau = sqrt(sigmaV2)))
  Kalman$set_observations(matrix(observations, ncol = 1))
  Kalman$first_step()
  for (t in 1:nobservations){
    Kalman$filtering_step(t-1)
  }
  loglikelihood = sum(Kalman$get_incremental_ll())
  if (log){
    return (model$dprior(rbind(theta,psi,sigmaV2),log = TRUE) + loglikelihood)
  } else {
    return (exp((model$dprior(rbind(theta,psi,sigmaV2),log = TRUE) + loglikelihood)))
  }
}
# Set parameters for Metropolis-Hastings
burnin = 3000
M = burnin + algorithmic_parameters$Ntheta
MH_cov = (cov.wt(t(thetas_smc[1:2,]),wt=smc_results$normw_history[[nobservations+1]])$cov)/5
thetas_MH = matrix(NA,nrow = 2,ncol = M)
thetas_MH[,1] = theta_star[1:2]
accepts = 0
# Track progress
print(paste("Started at:",Sys.time()))
progbar = txtProgressBar(min = 0,max = M,style=3)
time_start = proc.time()
# Run MH
for (i in 2:M){
  theta_new = fast_rmvnorm_transpose(1,thetas_MH[,i-1],MH_cov)
  if (model$dprior(theta_new) == -Inf){
    thetas_MH[,i] = thetas_MH[,i-1]
    setTxtProgressBar(progbar, i)
    next
  } else {
    logacceptance = dpost(theta_new) - dpost(thetas_MH[,i-1])
    logu = log(runif(1))
    if (logu <= logacceptance){
      accepts = accepts + 1
      thetas_MH[,i] = theta_new
    } else {
      thetas_MH[,i] = thetas_MH[,i-1]
    }
    setTxtProgressBar(progbar, i)
  }
}
# Track progress
time_end = proc.time()-time_start
print(time_end)
cat("acceptance rate = ", accepts/(M-1))
# Traceplots
par(mfrow=c(2,1))
index = (burnin+1):M
plot(index,thetas_MH[1,index],type='l')
plot(index,thetas_MH[2,index],type='l')
par(mfrow=c(1,1))
#--------------------------------------------------------------------------------------------
thetas_MH = thetas_MH[,index]

# Checking sample from the posterior distribution (marginal histogram)
Ntheta = algorithmic_parameters$Ntheta

g1 <- qplot(x = thetas_smc[1,], weight = normw_smc, geom = "blank") + geom_density(aes(y = ..density..))
g1 <- g1 + geom_density(data = data.frame(x = thetas_MH[1,]), aes(x = x, weight = NULL), col = "red")
g1 <- g1 + geom_density(data = data.frame(x = thetasPMMH[1,]), aes(x = x, weight = NULL), col = "blue")
g1

g2 <- qplot(x = thetas_smc[2,], weight = normw_smc, geom = "blank") + geom_density(aes(y = ..density..))
g2 <- g2 + geom_density(data = data.frame(x = thetas_MH[2,]), aes(x = x, weight = NULL), col = "red")
g2 <- g2 + geom_density(data = data.frame(x = thetasPMMH[2,]), aes(x = x, weight = NULL), col = "blue")
g2
