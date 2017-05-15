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
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$store_X = FALSE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$min_acceptance_rate = 0.45
algorithmic_parameters$nmoves = 2
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R

### Run SMC
kalman_module <<- Module( "kalman_mod", PACKAGE = "HyvarinenSSM")
smcKF_results <- smcKF_sampler(observations, model, algorithmic_parameters)
smc_results <- smc_sampler(observations, model, algorithmic_parameters)

### Run SMC_2 (non-tempered and tempered)
module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree
algorithmic_parameters$progress = TRUE
smc2_results <- hscore_continuous_no_tempering(observations, model, algorithmic_parameters)
algorithmic_parameters$progress = FALSE
model_withoutlikelihood = model
model_withoutlikelihood$likelihood = NULL
model_withoutlikelihood$dpredictive = NULL
smc2_results_temp <- hscore(observations, model_withoutlikelihood, algorithmic_parameters)

###############################################################################################
###############################################################################################
#### Sanity check: sample posterior via naive MH
# Computes the posterior density (target)
psi = model$psi
sigmaV2 = model$sigmaV2
dpost = function(theta, log = TRUE){
  Kalman <- new(kalman_module$Kalman)
  Kalman$set_parameters(list(rho = theta[1], sigma = sqrt(theta[2]), eta = psi, tau = sqrt(sigmaV2)))
  Kalman$set_observations(matrix(observations, ncol = 1))
  Kalman$first_step()
  for (t in 1:nobservations){
    Kalman$filtering_step(t-1)
  }
  loglikelihood = sum(Kalman$get_incremental_ll())
  if (log){
    return (model$dprior(theta,log = TRUE) + loglikelihood)
  } else {
    return (exp((model$dprior(theta,log = TRUE) + loglikelihood)))
  }
}
# Set parameters for Metropolis-Hastings
burnin = 500
M = burnin + algorithmic_parameters$Ntheta
MH_cov = (cov.wt(t(smcKF_results$thetas_history[[nobservations+1]]),wt=smcKF_results$normw_history[[nobservations+1]])$cov)/5
thetas_MH = matrix(NA,nrow = 2,ncol = M)
thetas_MH[,1] = c(0.7,1)
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

########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################
time_t = 50
#
thetas_MH = thetas_MH[,index]
#
thetas_smcKF <- smcKF_results$thetas_history[[time_t+1]]
normw_smcKF <- smcKF_results$normw_history[[time_t+1]]
#
thetas_smc <- smc_results$thetas_history[[time_t+1]]
normw_smc <- smc_results$normw_history[[time_t+1]]
#
thetas_smc2 <- smc2_results$thetas_history[[time_t+1]]
normw_smc2 <- smc2_results$weights_history[[time_t+1]]
#
thetas_smc2_temp <- smc2_results_temp$thetas_history[[time_t+1]]
normw_smc2_temp <- smc2_results_temp$normw_history[[time_t+1]]

# Checking sample from the posterior distribution (marginal histogram)
Ntheta = algorithmic_parameters$Ntheta
post = data.frame(from = factor(rep(c("smc","smcKF","smc2","smc2temp","MH"),each = Ntheta)))
post$theta1 = c(thetas_smc[1,],thetas_smcKF[1,],thetas_smc2[1,],thetas_smc2_temp[1,],thetas_MH[1,])
post$theta2 = c(thetas_smc[2,],thetas_smcKF[2,],thetas_smc2[2,],thetas_smc2_temp[2,],thetas_MH[2,])
post$weight = c(normw_smc,normw_smcKF,normw_smc2,normw_smc2_temp,rep(1/Ntheta,Ntheta))
#
plot_theta1 = ggplot(post) +  geom_density(aes(theta1, weight = weight, fill = from), alpha = 0.6) + theme(legend.position="none")
plot_theta2 = ggplot(post) +  geom_density(aes(theta2, weight = weight, fill = from), alpha = 0.6)
grid.arrange(plot_theta1, plot_theta2, ncol = 2, widths=c(1,1.35))

# Check the log-evidence (RESCALED BY 1/t)
results = data.frame(from = factor(rep(c("smc","smcKF","smc2","smc2temp"),each = time_t)))
results$time = rep(1:time_t, 4)
results$logevidence = c(smc_results$logevidence,smcKF_results$logevidence, smc2_results$logevidence, smc2_results_temp$logevidence)
ggplot(results) + geom_line(aes(time, logevidence/time, color = from), size = 1)

# # Checking sample from the posterior distribution (contour)
# ### Visual (qualitative) diagnostic
# M = 100
# x = seq(0.01, 1.5, length.out = M)
# y = seq(0.5,10, length.out = M)
# z = matrix(NA,ncol = M,nrow = M)
# print(paste("Started at:",Sys.time()))
# progbar = txtProgressBar(min = 0,max = M*M,style=3)
# count = 0
# time_start = proc.time()
# for (i in 1:M){
#   for (j in 1:M){
#     z[i,j] = dpost(c(x[i],y[j]))
#     count = count + 1
#     setTxtProgressBar(progbar, count)
#   }
# }
# time_end = proc.time()-time_start
# print(time_end)
# contour_df = expand.grid(x = x, y = y)
# contour_df$z = c(z)
# ggplot() +
#   geom_point(aes(thetas_smcKF[,1], thetas_smcKF[,2], alpha = normw_smcKF)) +
#   geom_contour(data = contour_df,aes(x,y,z=z,color = ..level..)) +
#   scale_colour_gradient(low="black", high="red") +
#   geom_point(aes(thetas_MH[index,1],thetas_MH[index,2]),color="yellow", size = 1, shape = 3) +
#   geom_point(aes(thetas_smc2[,1], thetas_smc2[,2], alpha = normw_smc2),color="purple", shape=15) +
#   geom_point(aes(x = theta_star[1], y = theta_star[2]), colour = "red", size = 10)
