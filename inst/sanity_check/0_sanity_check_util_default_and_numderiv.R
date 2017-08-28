##################################################################################################
# Sanity check for the function set_default_model (util_default.R).
# This checks that the automatic definitions of the missing fields of the model (e.g. using
# numerical differentiation) agree with the output using analytical definitions
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
nobservations = 15
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
algorithmic_parameters$Nx = 2^5
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$store_X = TRUE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$min_acceptance_rate = 0.2
algorithmic_parameters$nmoves = 2
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R
#--------------------------------------------------------------------------------------------
# Model with automatic predictive and numerical derivatives
model_likelihood_only = model
model_likelihood_only$dpredictive = NULL
model_likelihood_only$derivativelogdpredictive = NULL

# Model with automatic likelihood and numerical derivatives
model_predictive_only = model
model_predictive_only$likelihood = NULL
model_predictive_only$derivativelogdpredictive = NULL

# Model with no likelihood, no predictive, but exact derivatives
model_dobs = model
model_dobs$likelihood = NULL
model_dobs$dpredictive = NULL

#M odel with no likelihood, no predictive, and numerical derivatives
model_dobs_only = model
model_dobs_only$likelihood = NULL
model_dobs_only$dpredictive = NULL
model_dobs_only$derivativelogdpredictive = NULL
model_dobs_only$derivativelogdobs = NULL
#--------------------------------------------------------------------------------------------
model_all = list(model, model_likelihood_only, model_predictive_only, model_dobs, model_dobs_only)
results.df = data.frame()
post.df = data.frame()
for (i in 1:length(model_all)) {
  result = hscore(observations,model_all[[i]],algorithmic_parameters)
  results.df = rbind(results.df,data.frame(time = 1:nobservations,
                                           logevidence = result$logevidence,
                                           hscore = result$hscore,
                                           index = i))
  #########################################################################################
  ########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################
  #########################################################################################
  post.df = rbind(post.df,data.frame(theta1 = result$thetas_history[[nobservations+1]][1,],
                                     theta2 = result$thetas_history[[nobservations+1]][2,],
                                     theta3 = result$thetas_history[[nobservations+1]][3,],
                                     theta4 = result$thetas_history[[nobservations+1]][4,],
                                     weight = result$normw_history[[nobservations+1]],
                                     index = i))
}
#--------------------------------------------------------------------------------------------
# plot posterior marginals
plot_theta1 = ggplot(post.df) +  geom_density(aes(theta1, weight = weight, group = index, fill = factor(index)), alpha = 0.6) + theme(legend.position="none")
plot_theta2 = ggplot(post.df) +  geom_density(aes(theta2, weight = weight, group = index, fill = factor(index)), alpha = 0.6) + theme(legend.position="none")
plot_theta3 = ggplot(post.df) +  geom_density(aes(theta3, weight = weight, group = index, fill = factor(index)), alpha = 0.6) + theme(legend.position="none")
plot_theta4 = ggplot(post.df) +  geom_density(aes(theta4, weight = weight, group = index, fill = factor(index)), alpha = 0.6)
grid.arrange(plot_theta1, plot_theta2, plot_theta3, plot_theta4, ncol = 4, widths=c(1,1,1,1.5))
#--------------------------------------------------------------------------------------------
# Check the log-evidence (RESCALED BY 1/t)
ggplot(results.df) + geom_line(aes(time, logevidence/time, group = index, color = factor(index)), size = 1)
#--------------------------------------------------------------------------------------------
# Check the h-score (RESCALED BY 1/t)
ggplot(results.df) + geom_line(aes(time, hscore/time, group = index, color = factor(index)), size = 1)

