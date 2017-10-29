##################################################################################################
# Generic template to compute the hscore of a model
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(29)
#--------------------------------------------------------------------------------------------
#--------------------------------- TREE MODULE (for SMC2) -----------------------------------
#--------------------------------------------------------------------------------------------
# if SMC2 will be used, need tree module (NB: automatically called upon loading HyvarinenSSM package)
module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree
#--------------------------------------------------------------------------------------------
#---------------------------------        DATA            -----------------------------------
#--------------------------------------------------------------------------------------------
nobservations = 20
Y = rnorm(nobservations,0,1)
# observations in a matrix of dimensions dimy x nobservations
observations = matrix(Y, nrow = 1)
#--------------------------------------------------------------------------------------------
#---------------------------------        MODEL            ----------------------------------
#--------------------------------------------------------------------------------------------
model = list() # <<<<<<< DEFINE MODEL
# e.g. model = get_model_iid_gaussian_unknown_mean(0,1000)
#--------------------------------------------------------------------------------------------
#---------------------               ALGORITHMIC PARAMETERS         -------------------------
#--------------------------------------------------------------------------------------------
# Unspecified fields are set to their default values
# via set_default_algorithmic_parameters in util_default.R
algorithmic_parameters = list()
# algorithmic_parameters$Ntheta = ... # Number of particles theta
# algorithmic_parameters$ess_threshold = ... # ESS threshold for rejuvenation
# algorithmic_parameters$nmoves = ... # Number of moves per rejuvenation step
# algorithmic_parameters$Nx = ... # (ONLY for SMC2) Number of particles X
# algorithmic_parameters$adaptNx = ... # (ONLY for SMC2) Adaptive number of particles X
# algorithmic_parameters$min_acceptance_rate = ... # (ONLY if adaptNx) Acceptance rae threshold to increase Nx
# algorithmic_parameters$Nx_max = ... # (ONLY for SMC2 if adaptNx) Maximum value for Nx
##### complete list of algorithmic parameters can be found in util_default.R
#--------------------------------------------------------------------------------------------
#---------------------                HSCORE + LOGEVIDENCE          -------------------------
#--------------------------------------------------------------------------------------------
# Run SMC or SMC2
results = hscore(observations, model, algorithmic_parameters)
# exctract particles thetas and weights
if (results$algorithmic_parameters$store_theta) {
  thetas = results$thetas_history[[nobservations+1]]
  normw = results$normw_history[[nobservations+1]]
}
# exctract logevidence
logevidence = results$logevidence
# exctract hscore
if (results$algorithmic_parameters$hscore) {
  hscore = results$hscore
}
