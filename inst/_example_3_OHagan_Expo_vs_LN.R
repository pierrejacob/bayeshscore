##################################################################################################
# This implements example 3.3. in O'Hagan (1995).
# Model 1 = iid expo(theta) with reference prior
# Model 2 = iid Log-Normal(mu,sigma2) with reference prior
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
set.seed(1)
#--------------------------------------------------------------------------------------------
# Generate some data
nobservations = 300
Y = rexp(nobservations,1)
observations = matrix(Y, nrow = 1)# observations in a matrix of dimensions dimy by nobservations
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$ess_threshold = 0.5
#--------------------------------------------------------------------------------------------
# define models
# Note: set a and b close to 0, and sigma02 close to infinity for improper reference priors
model1 = get_model_iid_exponential(a = 1, b = 1)
model2 = get_model_iid_lognormal(mu0 = 0, sigma02 = 100, a = 1, b = 1)
models = list(model1, model2)
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------
results = data.frame()
for (m in 1:length(models)){
  result = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations, models[[m]], algorithmic_parameters)
  }
  for (r in 1:repl){
    results = rbind(results,data.frame(logevidence = result[[r]]$logevidence,
                                       hscore = result[[r]]$hscore,
                                       time = 1:nobservations,
                                       model = m,
                                       repl = r))
  }
}
#--------------------------------------------------------------------------------------------
# Check the log-evidence (WARNING: only meaningful for not so vague priors)
ggplot(results) +
  geom_line(aes(time, -logevidence/time, color = factor(model), group=interaction(model,repl)))
#--------------------------------------------------------------------------------------------
# Check the h-score
ggplot(results) +
  geom_line(aes(time, hscore/time, color = factor(model), group=interaction(model,repl)))
#--------------------------------------------------------------------------------------------
