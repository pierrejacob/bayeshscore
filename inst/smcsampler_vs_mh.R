##################################################################################################
# This checks that the outputs using SMC and SMC2 for continuous  observations
# match the exact results from Metropolis-Hastings in a linear gaussian example (with 2 parameters)
# where exact computation are achievable by using Kalman filters.
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(19)

# create data
nobservations = 10
model = get_model_lineargaussian()
theta_star = c(0.8,1,1,1)
sim = simulateData(model, theta = theta_star, nobservations)
X = sim$X
Y = sim$Y
observations = matrix(Y, nrow = model$dimY)
# observations in a matrix of dimensions dimY by nobservations


#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^14
algorithmic_parameters$Nx = 2^6
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$store_X = FALSE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$min_acceptance_rate = 0.5 # purposely set high to trigger some increase Nx steps
algorithmic_parameters$nmoves = 50
algorithmic_parameters$proposalmove <- get_independent_normal_proposal()
# algorithmic_parameters$proposalmove <- get_mixture_proposal()
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R

#--------------------------------------------------------------------------------------------
### Run SMC
smc_results = hscore(observations, model, algorithmic_parameters)

########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################
thetas_smc = smc_results$thetas_history[[nobservations+1]]
normw_smc = smc_results$normw_history[[nobservations+1]]
#--------------------------------------------------------------------------------------------
###########################################################################################
###########################################################################################
#### Sanity check: sample posterior via naive MH
# Computes the posterior density (target)
psi = model$psi
sigmaV2 = model$sigmaV2
kalman_module <<- Module( "kalman_mod", PACKAGE = "HyvarinenSSM")
dpost = function(theta, log = TRUE){
  Kalman <<- new(kalman_module$Kalman)
  Kalman$set_parameters(list(rho = theta[1], sigma = sqrt(theta[2]), eta = theta[3], tau = sqrt(theta[4])))
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
# M = burnin + algorithmic_parameters$Ntheta
M = burnin + 5000
MH_cov = (cov.wt(t(thetas_smc),wt=smc_results$normw_history[[nobservations+1]])$cov)/5
thetas_MH = matrix(NA,nrow = 4,ncol = M)
thetas_MH[,1] = theta_star
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
par(mfrow=c(2,2))
index = (burnin+1):M
plot(index,thetas_MH[1,index],type='l')
plot(index,thetas_MH[2,index],type='l')
plot(index,thetas_MH[3,index],type='l')
plot(index,thetas_MH[4,index],type='l')
par(mfrow=c(1,1))
#--------------------------------------------------------------------------------------------
thetas_MH = thetas_MH[,index]

# Checking sample from the posterior distribution (marginal histogram)
Ntheta = algorithmic_parameters$Ntheta

dim(thetas_smc)

g1 <- qplot(x = thetas_smc[1,], weight = normw_smc, geom = "blank") + geom_density(aes(y = ..density..))
g1 <- g1 + geom_density(data = data.frame(x = thetas_MH[1,]), aes(x = x, weight = NULL), col = "red")
g1


g2 <- qplot(x = thetas_smc[2,], weight = normw_smc, geom = "blank") + geom_density(aes(y = ..density..), adjust = 0.1)
g2 <- g2 + geom_density(data = data.frame(x = thetas_MH[2,]), aes(x = x, weight = NULL), col = "red", adjust = 0.1)
g2

g2 <- qplot(x = thetas_smc[2,], weight = normw_smc, geom = "blank") + geom_histogram(aes(y = ..density..), alpha = 0.5, binwidth = 0.1)
g2 <- g2 + geom_histogram(data = data.frame(x = thetas_MH[2,]), aes(x = x, weight = NULL, y = ..density..), fill = "red", alpha = 0.5, binwidth = 0.1)
g2


g3 <- qplot(x = thetas_smc[3,], weight = normw_smc, geom = "blank") + geom_histogram(aes(y = ..density..), alpha = 0.5, binwidth = 0.05)
g3 <- g3 + geom_histogram(data = data.frame(x = thetas_MH[3,]), aes(x = x, weight = NULL, y = ..density..), fill = "red", alpha = 0.5, binwidth = 0.05)
g3

g4 <- qplot(x = thetas_smc[4,], weight = normw_smc, geom = "blank") + geom_histogram(aes(y = ..density..), alpha = 0.5, binwidth = 0.05)
g4 <- g4 + geom_histogram(data = data.frame(x = thetas_MH[4,]), aes(x = x, weight = NULL, y = ..density..), fill = "red", alpha = 0.5, binwidth = 0.05)
g4


# post = data.frame(from = factor(rep(c("smc", "MH"),each = Ntheta)))
# post$theta1 = c(thetas_smc[1,],thetas_MH[1,])
# post$theta2 = c(thetas_smc[2,],thetas_MH[2,])
# post$weight = c(normw_smc,rep(1/Ntheta,Ntheta))
# # plot posterior marginals
# plot_theta1 = ggplot(post) +  geom_density(aes(theta1, weight = weight, fill = from), alpha = 0.6) + theme(legend.position="none")
# plot_theta2 = ggplot(post) +  geom_density(aes(theta2, weight = weight, fill = from), alpha = 0.6)
# grid.arrange(plot_theta1, plot_theta2, ncol = 2, widths=c(1,1.35))
#
# # #--------------------------------------------------------------------------------------------
# # # Check the log-evidence (RESCALED BY 1/t)
# # results = data.frame(from = factor(rep(c("smc","smc2"),each = nobservations)))
# # results$time = rep(1:nobservations, 2)
# # results$logevidence = c(smc_results$logevidence, smc2_results$logevidence)
# # ggplot(results) + geom_line(aes(time, logevidence/time, color = from), size = 1)
# #
# # #--------------------------------------------------------------------------------------------
# # # Check the hscore (RESCALED BY 1/t)
# # results$hscore = c(smc_results$hscore, smc2_results$hscore)
# # ggplot(results) + geom_line(aes(time, hscore/time, color = from), size = 1)
