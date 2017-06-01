##################################################################################################
# This checks that the outputs from SMC and SMC2 agree in a linear gaussian case (4 parameters)
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(19)
#--------------------------------------------------------------------------------------------
model = get_model_lineargaussian()
#--------------------------------------------------------------------------------------------
# create data
nobservations = 15
theta_star = c(0.8,1,0.5,1)
sim = simulateData(model, theta = theta_star, nobservations)
X = sim$X
Y = sim$Y
observations = matrix(Y, nrow = model$dimY)
# observations in a matrix of dimensions dimY by nobservations
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^12
algorithmic_parameters$Nx = 2^5
# algorithmic_parameters$Nx_max = 2^6
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$store_X = TRUE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$min_acceptance_rate = 0.5 # purposely set high to trigger increase Nx step and test Nx_max
algorithmic_parameters$nmoves = 5 # purposely set high for sanity check
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R
#--------------------------------------------------------------------------------------------
### Run SMC
smc_results = hscore(observations, model, algorithmic_parameters)
### Run SMC_2
model_withoutlikelihood = model
model_withoutlikelihood$likelihood = NULL # this forces the use of SMC2
model_withoutlikelihood$dpredictive = NULL # this forces the use of SMC2
smc2_results = hscore(observations, model_withoutlikelihood, algorithmic_parameters)
#--------------------------------------------------------------------------------------------
#########################################################################################
########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################
thetas_smc = smc_results$thetas_history[[nobservations+1]]
normw_smc = smc_results$normw_history[[nobservations+1]]
#
thetas_smc2 = smc2_results$thetas_history[[nobservations+1]]
normw_smc2 = smc2_results$normw_history[[nobservations+1]]
#--------------------------------------------------------------------------------------------
# Checking sample from the posterior distribution (marginal histogram)
Ntheta = algorithmic_parameters$Ntheta
post = data.frame(from = factor(rep(c("smc","smc2"),each = Ntheta)))
post$theta1 = c(thetas_smc[1,],thetas_smc2[1,])
post$theta2 = c(thetas_smc[2,],thetas_smc2[2,])
post$theta3 = c(thetas_smc[3,],thetas_smc2[3,])
post$theta4 = c(thetas_smc[4,],thetas_smc2[4,])
post$weight = c(normw_smc,normw_smc2)
# plot posterior marginals
plot_theta1 = ggplot(post) +  geom_density(aes(theta1, weight = weight, fill = from), alpha = 0.6) + theme(legend.position="none")
plot_theta2 = ggplot(post) +  geom_density(aes(theta2, weight = weight, fill = from), alpha = 0.6) + theme(legend.position="none")
plot_theta3 = ggplot(post) +  geom_density(aes(theta3, weight = weight, fill = from), alpha = 0.6) + theme(legend.position="none")
plot_theta4 = ggplot(post) +  geom_density(aes(theta4, weight = weight, fill = from), alpha = 0.6)
grid.arrange(plot_theta1, plot_theta2, plot_theta3, plot_theta4, ncol = 4, widths=c(1,1,1,1.5))
#--------------------------------------------------------------------------------------------
# Check the log-evidence (RESCALED BY 1/t)
results = data.frame(from = factor(rep(c("smc","smc2"),each = nobservations)))
results$time = rep(1:nobservations, 2)
results$logevidence = c(smc_results$logevidence,smc2_results$logevidence)
ggplot(results) + geom_line(aes(time, logevidence/time, color = from), size = 1) +
  geom_vline(data=data.frame(rejuvenationtime = smc2_results$rejuvenation_times,
                             i = 1:length(smc2_results$rejuvenation_times),
                             incrNx = factor(smc2_results$rejuvenation_times %in% smc2_results$increase_Nx_times)),
             aes(xintercept = rejuvenationtime, group = i, linetype = incrNx))
#--------------------------------------------------------------------------------------------
# Check the h-score (RESCALED BY 1/t)
results$hscore = c(smc_results$hscore,smc2_results$hscore)
ggplot(results) + geom_line(aes(time, hscore/time, color = from), size = 1) +
  geom_vline(data=data.frame(rejuvenationtime = smc2_results$rejuvenation_times,
                             i = 1:length(smc2_results$rejuvenation_times),
                             incrNx = factor(smc2_results$rejuvenation_times %in% smc2_results$increase_Nx_times)),
             aes(xintercept = rejuvenationtime, group = i, linetype = incrNx))


#--------------------------------------------------------------------------------------------
# estimate (marginal) filtering means from SMC2
xmean = rep(0,nobservations)
for (t in 1:nobservations) {
  PF = smc2_results$PF_history[[t+1]]
  normw = smc2_results$normw_history[[t+1]]
  xmean[t] = sum(sapply(1:Ntheta,function(i)sum(PF[[i]]$X*PF[[i]]$xnormW)*normw[i]))
}
# estimate (marginal) filtering means via KF (using the same particles thetas from SMC2)
xmeanKF = rep(0,nobservations)
for (t in 1:nobservations) {
  normw = smc2_results$normw_history[[t+1]]
  thetas = smc2_results$thetas_history[[t+1]]
  for (i in 1:Ntheta){
    phi = thetas[1,i]
    sigmaW2 = thetas[2,i]
    psi = thetas[3,i]
    sigmaV2 = thetas[4,i]
    initial_mean = 0
    initial_var = (sigmaW2)/(1-phi^2)
    KF = KF_filtering(observations[,1:t,drop=FALSE],phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var)
    xmeanKF[t] = xmeanKF[t] + KF[[t]]$muX_t_t*normw[i]
  }
}
# plot filtering means
ggplot(data.frame(time = rep(1:nobservations,2), marginalfiltermean = c(xmean,xmeanKF),
                  from = factor(rep(c("smc2","KF"),each=nobservations)))) +
  geom_line(aes(time, marginalfiltermean, color = from),size=1) +
  geom_vline(data=data.frame(rejuvenationtime = smc2_results$rejuvenation_times,
                             i = 1:length(smc2_results$rejuvenation_times),
                             incrNx = factor(smc2_results$rejuvenation_times %in% smc2_results$increase_Nx_times)),
             aes(xintercept = rejuvenationtime, group = i, linetype = incrNx))


# result = data.frame()
# for (t in 1:nobservations) {
#   PF = smc2_results$PF_history[[t+1]]
#   for (i in 1:Ntheta){
#     result = rbind(result, data.frame(X = c(PF[[i]]$X), WX = c(PF[[i]]$xnormW), index = i, time = t))
#   }
# }
#
# # Rejuvenation at step 6 with no Nx increase
# time_min = min(nobservations,5)
# time_max = min(nobservations,7)
# ggplot(subset(result,time_min<=time & time<=time_max)) +
#   geom_density(aes(X,weight=WX,group=interaction(index,time),fill=time))
#
#
# # Rejuvenation at step 13 with Nx increase
# time_min = min(nobservations,12)
# time_max = min(nobservations,14)
# ggplot(subset(result,time_min<=time & time<=time_max)) +
#   geom_density(aes(X,weight=WX,group=interaction(index,time),fill=time))
