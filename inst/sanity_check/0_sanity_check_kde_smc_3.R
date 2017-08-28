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
nu0 = 1
s02 = 1
nobservations = 10
observations = matrix(rnorm(nobservations), nrow = 1)# observations in a matrix of dimensions dimy x nobservations
model = get_model_iid_gaussian_unknown_variance(nu0,s02)
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_thetas_history = TRUE

algorithmic_parameters$use_kde = TRUE
algorithmic_parameters$kde_opt = list(Ny = 10^6, nb_steps = 2)
#--------------------------------------------------------------------------------------------
### Run SMC
results = hscore(observations, model, algorithmic_parameters)
#--------------------------------------------------------------------------------------------
#########################################################################################
#########################################################################################
########### BE CAREFUL, SMC starts with the prior sample at t = 1 #######################
#----------------------------------------------------------------------------------------
hscore_exact = rep(NA,nobservations)
hscore_numderiv = rep(NA,nobservations)

for (t in 1:nobservations){
  nu_t = nu0 + (t-1)
  if (t == 1){st2 = nu0/nu_t}
  if (t >= 2){st2 = (nu0+sum(observations[,1:(t-1)]^2))/nu_t}
  # numDeriv log-derivatives (as a sanity check of the previous analytical expressions)
  ll_exact_fn = function(x) {dtscaled(x,nu_t,st2,TRUE)}
  # compute h score exact
  s = sum(observations[,1:t]^2)
  hscore_exact[t] = ((nu0+t)/((nu0+s)^2))*((nu0+t+4)*observations[,t]^2-2*(nu0+s))
  # compute h score exact via numDeriv
  d1ll_t = grad(ll_exact_fn,observations[,t])
  d2ll_t = hessian(ll_exact_fn,observations[,t])
  hscore_numderiv[t] = (2*d2ll_t + (d1ll_t)^2)
}


hscoreKDE = results$hscoreKDE
if (algorithmic_parameters$kde_opt$nb_steps < nobservations) {
  lastindex = algorithmic_parameters$kde_opt$nb_steps
  hscoreKDE = c(hscoreKDE[1:lastindex], hscoreKDE[lastindex] + cumsum(diff(results$hscore)[lastindex:(nobservations-1)]))
}


ggplot() +
  geom_line(aes(1:nobservations, cumsum(hscore_numderiv)/(1:nobservations)),color="red") +
  geom_line(aes(1:nobservations, hscoreKDE/(1:nobservations)),color="blue",size = 1) +
  geom_line(aes(1:nobservations, cumsum(hscore_exact)/(1:nobservations)),color="red",size=2,linetype=2)





diff(results$hscore)[algorithmic_parameters$kde_opt$nb_steps:(nobservations-1)]


