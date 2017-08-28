##################################################################################################
# This checks that the results are saved properly and can be used to resume a run
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(19)
module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree
#--------------------------------------------------------------------------------------------
# create data
nobservations = 100
model = get_model_lineargaussian()
theta_star = c(0.8,1,1,1)
#--------------------------------------------------------------------------------------------
sim = simulateData(model, theta = theta_star, nobservations)
X = sim$X
Y = sim$Y
observations = matrix(Y, nrow = model$dimY)
# observations in a matrix of dimensions dimY by nobservations
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$Nx = 2^10
algorithmic_parameters$Nx_max = 2^10
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_thetas_history = TRUE
algorithmic_parameters$store_X_history = TRUE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$min_acceptance_rate = 0.2 # purposely set high to trigger increase Nx step and test Nx_max
algorithmic_parameters$nmoves = 2
algorithmic_parameters$save = TRUE
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R
#--------------------------------------------------------------------------------------------
#########################################################################################
#########################################################################################
#########################################################################################
#############               PURPOSELY INTERRUPT COMPUTATION               ###############
#########################################################################################
#########################################################################################
#########################################################################################
#
# WARNING: the save must be an RDS file (extension .rds)
#
### Run SMC0
setTimeLimit(elapsed = 20) # purposely set a time budget
algorithmic_parameters$savefilename = "partial_smc_results.rds"
tryCatch(hscore(observations,model,algorithmic_parameters), error = function(e) cat("--- Time out ---\n"))
### Run SMC_2
model_nolikelihood = model
model_nolikelihood$likelihood = NULL # this forces the use of SMC2
model_nolikelihood$dpredictive = NULL # this forces the use of SMC2
setTimeLimit(elapsed = 60) # purposely set a time budget
algorithmic_parameters$savefilename = "partial_smc2_results.rds"
tryCatch(hscore(observations,model_nolikelihood,algorithmic_parameters), error = function(e) cat("--- Time out ---\n"))
#--------------------------------------------------------------------------------------------
setTimeLimit() # This resets the time limit to infinity
#########################################################################################
########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################
#########################################################################################
partial_smc_results = readRDS("partial_smc_results.rds")
partial_smc2_results = readRDS("partial_smc2_results.rds")
#
nobs_smc = partial_smc_results$t
nobs_smc2 = partial_smc2_results$t
partial_nobs = min(nobs_smc, nobs_smc2)
#
thetas_smc = partial_smc_results$thetas_history[[partial_nobs+1]]
normw_smc = partial_smc_results$normw_history[[partial_nobs+1]]
#
thetas_smc2 = partial_smc2_results$thetas_history[[partial_nobs+1]]
normw_smc2 = partial_smc2_results$normw_history[[partial_nobs+1]]
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
results = data.frame(from = factor(rep(c("smc","smc2"),each = partial_nobs)))
results$time = rep(1:partial_nobs, 2)
results$logevidence = c(cumsum(partial_smc_results$incr_logevidence[1:partial_nobs]),
                        cumsum(partial_smc2_results$incr_logevidence[1:partial_nobs]))
ggplot(results) + geom_line(aes(time, logevidence/time, color = from), size = 1)
#--------------------------------------------------------------------------------------------
# Check the h-score (RESCALED BY 1/t)
results$hscore = c(cumsum(partial_smc_results$incr_hscore[1:partial_nobs]),
                        cumsum(partial_smc2_results$incr_hscore[1:partial_nobs]))
ggplot(results) + geom_line(aes(time, hscore/time, color = from), size = 1)


#--------------------------------------------------------------------------------------------
# estimate (marginal) filtering means from SMC2
xmean = rep(0,nobs_smc2)
Ntheta = partial_smc2_results$algorithmic_parameters$Ntheta
for (t in 1:nobs_smc2) {
  PF = partial_smc2_results$PF_history[[t+1]]
  normw = partial_smc2_results$normw_history[[t+1]]
  xmean[t] = sum(sapply(1:Ntheta,function(i)sum(PF[[i]]$X*PF[[i]]$xnormW)*normw[i]))
}
# estimate (marginal) filtering means via KF (using the same particles thetas from SMC2)
xmeanKF = rep(0,nobs_smc2)
for (t in 1:nobs_smc2) {
  normw = partial_smc2_results$normw_history[[t+1]]
  thetas = partial_smc2_results$thetas_history[[t+1]]
  for (i in 1:Ntheta){
    phi = thetas[1,i]
    sigmaW2 = thetas[2,i]
    psi = thetas[3,i]
    sigmaV2 = thetas[4,i]
    initial_mean = 0
    initial_var = (sigmaW2)/(1-phi^2)
    KF = KF_filtering(partial_smc2_results$observations[,1:t,drop=FALSE],phi,psi,sigmaV2,sigmaW2,initial_mean,initial_var)
    xmeanKF[t] = xmeanKF[t] + KF[[t]]$muX_t_t*normw[i]
  }
}
# plot filtering means
ggplot(data.frame(time = rep(1:nobs_smc2,2), marginalfiltermean = c(xmean,xmeanKF),
                  from = factor(rep(c("smc2","KF"),each=nobs_smc2)))) +
  geom_line(aes(time, marginalfiltermean, color = from),size=1)

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#-----------------------  Check that the trees were properly saved --------------------------
#--------------------------------------------------------------------------------------------
# check the content of the RDS save file
print(names(partial_smc2_results))
# Reconstruct PF_history
PF_history = lapply(1:partial_smc2_results$t,
                    function(t)lapply(1:length(partial_smc2_results$trees_attributes_history[[t+1]]),
                                      function(i)c(partial_smc2_results$PF_history_no_tree[[t+1]][[i]],
                                                   tree = tree_reconstruct(partial_smc2_results$trees_attributes_history[[t+1]][[i]]))))
PF_history = c(list(NULL),PF_history)
# Plot all paths from one of the latest trees (for some arbitrary theta)
path.df = data.frame(x = c(sapply(1:PF_history[[nobs_smc2+1]][[1]]$Nx,function(i)PF_history[[nobs_smc2+1]][[1]]$tree$get_path(i-1))),
                     time = rep(1:nobs_smc2,PF_history[[nobs_smc2+1]][[1]]$Nx),
                     index = rep(1:PF_history[[nobs_smc2+1]][[1]]$Nx,each=nobs_smc2))
ggplot(path.df) + geom_line(aes(time,x,group=index))


#----------------------------------------------------------------------------------------------------------------#
# # NOTE: lapply vs for loop to reconstruct trees
# lapply_version = function() {
#   PF_history = lapply(1:partial_smc2_results$t,
#                       function(t)lapply(1:length(partial_smc2_results$trees_attributes_history[[t+1]]),
#                                         function(i)c(partial_smc2_results$PF_history_no_tree[[t+1]][[i]],
#                                                      tree = tree_reconstruct(partial_smc2_results$trees_attributes_history[[t+1]][[i]]))))
# }
# for_version = function() {
#   PF_history = vector("list",partial_smc2_results$t)
#   PFs = vector("list",Ntheta)
#   for (t in 1:partial_smc2_results$t){
#     for (i in 1:Ntheta){
#       PFs[[i]] = c(partial_smc2_results$PF_history_no_tree[[t+1]][[i]],
#                    tree = tree_reconstruct(partial_smc2_results$trees_attributes_history[[t+1]][[i]]))
#     }
#     PF_history[[t]] = PFs
#   }
# }
# microbenchmark::microbenchmark(lapply_version(),for_version(),times = 20)
