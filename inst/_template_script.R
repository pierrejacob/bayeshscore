##################################################################################################
# Generic template to compute the hscore of a model
##################################################################################################
library(bayeshscore)
set.seed(29)
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
#>>>>>>>>> complete description of fields can be found in _template_model.R
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
# algorithmic_parameters$adaptNx = ... # (ONLY for SMC2) Allow adaptive number of particles X
# algorithmic_parameters$min_acceptance_rate = ... # (ONLY if adaptNx) Acceptance rae threshold to increase Nx
# algorithmic_parameters$Nx_max = ... # (ONLY if adaptNx) Maximum value for Nx
#>>>>>>>>> complete description of algorithmic parameters can be found in util_default.R
#--------------------------------------------------------------------------------------------
#---------------------                HSCORE + LOGEVIDENCE          -------------------------
#--------------------------------------------------------------------------------------------
# Compute the hscore (the function hscore is a wrapper that either calls smc or smc2 depending on whether
# the likelihood is available)
# NB: the log-evidence is also computed on the fly
results = hscore(observations, model, algorithmic_parameters)
#--------------------------------------------------------------------------------------------
# # an equivalent result would be obtained by the following calls
# model = set_default_model(model)
# algorithmic_parameters = set_default_algorithmic_parameters(observations,model,algorithmic_parameters)
# algorithmic_parameters$hscore = TRUE
# # run whichever line is appropriate (smc for tractable likelihood, smc2 otherwise)
# results = smc(observations, model, algorithmic_parameters)
# results = smc2(observations, model, algorithmic_parameters)
#--------------------------------------------------------------------------------------------
# results is then a list made of the following objects (some are set to NULL depending on what the user asks for)
#>>> thetas = last set of particles thetas (dim_theta by Ntheta matrix)
#>>> normw = associated normalized weights (vector of length Ntheta)
#>>> byproducts or PFs = last list of byproducts (e.g. particle filters in the case of smc2)
#>>> logtargetdensities = last target log-densities
#>>> thetas_history = list of all successive sets of particles thetas (time 1 = prior, so length = nobservations + 1)
#>>> normw_history = list of associated normalized weights
#>>> logtargetdensities_history = list of associated target log-densities
#>>> byproducts_history of PF_history = list of associated byproducts (e.g. particle filters in the case of smc2)
#>>> logevidence = cumulative logevidence
#>>> hscore = cumulative hscore (using Fisher/Louis type identities)
#>>> hscoreDDE = cumulative hscore (using kernel density estimation)
#>>> ESS = successive ESS
#>>> rejuvenation_times = successive times when rejuvenation occured
#>>> rejuvenation_rate = associated acceptance rates
#>>> method = method used ('SMC' or 'SMC2')
#>>> algorithmic_parameters = list of algorithmic parameters used

