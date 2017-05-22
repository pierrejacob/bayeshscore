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
model = get_model_iid_gaussian_unknown_mean(0,1000) # <<<<<<< DEFINE MODEL
#--------------------------------------------------------------------------------------------
#---------------------               ALGORITHMIC PARAMETERS         -------------------------
#--------------------------------------------------------------------------------------------
# Unspecified fields are set to their default values (in parentheses)
# via set_default_algorithmic_parameters in util_default.R
algorithmic_parameters = list()
# algorithmic_parameters$Ntheta = 2^10 # <<<< Number of particles theta (2^7)
# algorithmic_parameters$ess_threshold = ... # <<<< ESS threshold for rejuvenation (0.5)
# algorithmic_parameters$nmoves = ... # <<<< Number of moves per rejuvenation step (1)
# algorithmic_parameters$Nx = ... # <<<< (ONLY for SMC2) Number of particles X (big O of nobservations)
# algorithmic_parameters$adaptNx = ... # <<<< (ONLY for SMC2) Adaptive number of particles X (TRUE)
# algorithmic_parameters$min_acceptance_rate = ... # <<<< (ONLY for SMC2 if adaptNx) Acceptance rae threshold to increase Nx (0.2)
# algorithmic_parameters$Nx_max = ... # <<<< (ONLY for SMC2 if adaptNx) Maximum value for Nx
# algorithmic_parameters$hscore = ... # <<<< if TRUE, computes the prequential hscore (TRUE)
# algorithmic_parameters$verbose = ... # <<<< if TRUE, displays tempering steps and acceptance rates (TRUE)
# algorithmic_parameters$progress = ... # <<<< if TRUE, displays progress bar (FALSE)
# algorithmic_parameters$store_theta = ... # <<<< if TRUE, stores the history of particles thetas (TRUE for smc and smc2, FALSE for hscore)
# algorithmic_parameters$store_X = ... # <<<< (ONLY for SMC2) if TRUE, stores the history of particles X (FALSE)
# algorithmic_parameters$store_byproducts = ... # <<<< (ONLY for SMC) if TRUE, stores the history of byproducts (FALSE)
# algorithmic_parameters$save = ... # <<<< if TRUE, saves intermediary results in savefilename (FALSE)
# algorithmic_parameters$savefilename = ... # <<<< (ONLY if save) filename to save intermediary results
# algorithmic_parameters$time_budget = ... # <<<< time budget for computation in number of seconds (NULL)
# algorithmic_parameters$resampling = ... # <<<< resampling function (systematic resampling)
# algorithmic_parameters$proposalmove = ... # <<<< proposal for rejuvenation moves (get_proposal_mixture())
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
