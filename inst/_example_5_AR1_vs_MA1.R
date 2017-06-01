##################################################################################################
# Example 5: AR(1) vs MA(1)
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
set.seed(19)

# Define model
nu0 = 1
sigma02 = 1
nb_models = 2
model = function(i){
  if (i==1){return(get_model_ARMA(1,0,nu0, sigma02))} # AR(1)
  if (i==2){return(get_model_ARMA(0,1,nu0, sigma02))} # MA(1)
}

# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$verbose = TRUE
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------


##################################################################################################
# Case 1: true model = AR(1)
##################################################################################################
nobservations = 20
true_model = 1
true_theta = c(0.5,2)
observations = simulateData(model(true_model),true_theta,nobservations)$Y
# observations in a matrix of dimensions dimy x nobservations
#--------------------------------------------------------------------------------------------
results1_all = data.frame()
### Compute logevidence and hscore for each model
for (m in 1:nb_models){
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations, model(m), algorithmic_parameters)
  }
  for (r in 1:repl){
    results1_all = rbind(results1_all,data.frame(logevidence = results[[r]]$logevidence,
                                                 hscore = results[[r]]$hscore,
                                                 time = 1:nobservations,
                                                 model = factor(m),
                                                 repl = r))
  }
}
### Plot log evidence
ggplot(results1_all) +
  geom_line(aes(time,-logevidence/time,color=model,group=interaction(model,repl)))
### Plot H score
ggplot(results1_all) +
  geom_line(aes(time,hscore/time,color=model,group=interaction(model,repl)))

##################################################################################################
# Case 2: true model = MA(1)
##################################################################################################
nobservations = 20
true_model = 2
true_theta = c(0.5,2)
observations = simulateData(model(true_model),true_theta,nobservations)$Y
# observations in a matrix of dimensions dimy x nobservations
#--------------------------------------------------------------------------------------------
results2_all = data.frame()
### Compute logevidence and hscore for each model
for (m in 1:nb_models){
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations, model(m), algorithmic_parameters)
  }
  for (r in 1:repl){
    results2_all = rbind(results2_all,data.frame(logevidence = results[[r]]$logevidence,
                                                 hscore = results[[r]]$hscore,
                                                 time = 1:nobservations,
                                                 model = factor(m),
                                                 repl = r))
  }
}
### Plot log evidence
ggplot(results2_all) +
  geom_line(aes(time,-logevidence/time,color=model,group=interaction(model,repl)))
### Plot H score
ggplot(results2_all) +
  geom_line(aes(time,hscore/time,color=model,group=interaction(model,repl)))


##################################################################################################
# Compare model 1 and model 2 in each case
##################################################################################################

