##################################################################################################
# This checks that the outputs using SMC and SMC2 match the exact results
# in an iid Normal case (with conjugate prior so that everything can be computed analytically)
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
set.seed(29)

# Define model and data
nobservations = 5
observations = matrix(rnorm(nobservations), nrow = 1)# observations in a matrix of dimensions dimy x nobservations
model = get_model_SVLevy_singlefactor(1:nobservations)
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_thetas_history = TRUE

algorithmic_parameters$use_kde = TRUE
algorithmic_parameters$kde_opt = list(Ny = 10^4, nb_steps = nobservations)
#--------------------------------------------------------------------------------------------
### Run SMC2
results = hscore(observations, model, algorithmic_parameters)
#--------------------------------------------------------------------------------------------
#########################################################################################
#########################################################################################
########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################
#----------------------------------------------------------------------------------------


hscoreKDE = results$hscoreKDE
if (algorithmic_parameters$kde_opt$nb_steps < nobservations) {
  lastindex = algorithmic_parameters$kde_opt$nb_steps
  hscoreKDE = c(hscoreKDE[1:lastindex], hscoreKDE[lastindex] + cumsum(diff(results$hscore)[lastindex:(nobservations-1)]))
}

ggplot() +
  geom_line(aes(1:nobservations, results$hscore/(1:nobservations)),color="red") +
  geom_line(aes(1:nobservations, hscoreKDE/(1:nobservations)),color="blue",size = 1)

