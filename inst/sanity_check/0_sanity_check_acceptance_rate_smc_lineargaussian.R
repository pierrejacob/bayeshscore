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
algorithmic_parameters$nmoves = 10 # purposely set high to check the coherence of acceptance rates
algorithmic_parameters$proposalmove = get_proposal_mixture()
#--------------------------------------------------------------------------------------------
# Run SMC
smc_result = hscore(observations, model, algorithmic_parameters)
#--------------------------------------------------------------------------------------------
