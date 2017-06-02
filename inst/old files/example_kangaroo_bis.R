rm(list = ls())
library(doParallel)
library(HyvarinenSSM)
library(gridExtra)
library(wesanderson)
set.seed(19)

#=======================================================================
#=======================================================================
#=======================================================================
#====== WARNING: current crashes when trying to use parallel foreach ===
#=======================================================================
#========= (not sure how to properly export the Cpp Class TreeClass ) ==
#=======================================================================


# Define model and data
model1 = get_model_kangarooLogistic()
model2 = get_model_kangarooExponential()
model3 = get_model_kangarooRandomwalk()
# all_models = list(model1,model2,model3)
all_models = list(model2,model3)
dataset = data_kangaroo


nobservations = 3
observations = dataset[1:2,1:nobservations]
# nobservations = ncol(observations)


# set algorithmic parameters
algorithmic_parameters <- list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$Nx = 2^7
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$store_X = FALSE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$min_acceptance_rate = 0.1
algorithmic_parameters$nmoves = 2
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R

module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree

M = all_models[[1]]
M$observation_type = 'discrete'
results = hscore(observations, M, algorithmic_parameters)
