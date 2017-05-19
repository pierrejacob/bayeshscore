##################################################################################################
# This checks the acceptance rates of the SMC on a linear gaussian case
# For this case: the SMC uses the byproduct feature (with Kalman filters)
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
set.seed(29)
module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree
#--------------------------------------------------------------------------------------------
model = get_model_lineargaussian()
#--------------------------------------------------------------------------------------------
# create data
nobservations = 20
theta_star = c(0.8,1,1,1)
sim = simulateData(model, theta = theta_star, nobservations)
X = sim$X
Y = sim$Y
observations = matrix(Y, nrow = model$dimY)
# observations in a matrix of dimensions dimY by nobservations
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$nmoves = 20 # purposely set high to check the coherence of acceptance rates
algorithmic_parameters$proposalmove = get_proposal_mixture()
#--------------------------------------------------------------------------------------------
# Remove likelihood and predictive in order to force the use of SMC2
model_nolikelihood = model
model_nolikelihood$likelihood = NULL
model_nolikelihood$dpredictive = NULL
# Run SMC2
smc2_result = hscore(observations, model_nolikelihood, algorithmic_parameters)
#--------------------------------------------------------------------------------------------
