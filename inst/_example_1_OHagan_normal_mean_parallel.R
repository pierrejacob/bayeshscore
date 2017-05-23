##################################################################################################
# This implements example 3.1. in O'Hagan (1995).
# The model is N(theta,1) with conjugate prior on theta, i.e. theta follows N(0,sigma2prior).
# We compute the logevidence and the prequential Hyvarinen score for increasing vagueness
# (i.e. increasing values of sigma2prior).
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
set.seed(29)

# Define model and data
nobservations = 20
Y = rnorm(nobservations,0,1)
observations = matrix(Y, nrow = 1)# observations in a matrix of dimensions dimy x nobservations

#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
algorithmic_parameters$ess_threshold = 0.5
algorithmic_parameters$nmoves = 2
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------
sigma2prior_all = c(100,1000,10000,100000)
results_all = data.frame()
post_all = data.frame()
#--------------------------------------------------------------------------------------------
### Compute logevidence and hscore for increasing vagueness
for (s in 1:length(sigma2prior_all)){
  sigma2prior = sigma2prior_all[s]
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations, get_model_iid_gaussian_unknown_mean(0,sigma2prior), algorithmic_parameters)
  }
  for (r in 1:repl){
    results_all = rbind(results_all,data.frame(logevidence = results[[r]]$logevidence,
                                               hscore = results[[r]]$hscore,
                                               time = 1:nobservations,
                                               sigma2prior = sigma2prior,
                                               repl = r))
    post_all = rbind(post_all,data.frame(theta = c(results[[r]]$thetas_history[[nobservations+1]]),
                                                   W = results[[r]]$normw_history[[nobservations+1]],
                                                   sigma2prior = sigma2prior,
                                                   repl = r))
  }
}
#--------------------------------------------------------------------------------------------
# Checking sample from the posterior distribution (marginal histogram)
ggplot(post_all) +
  geom_density(aes(theta,weight=W,fill=factor(sigma2prior),group=interaction(sigma2prior,repl)),alpha=0.6)
#--------------------------------------------------------------------------------------------
# Check the log-evidence
ggplot(results_all) +
  geom_line(aes(time, -logevidence, color = factor(sigma2prior),group=interaction(sigma2prior,repl)))
#--------------------------------------------------------------------------------------------
# Check the h-score
ggplot(results_all) +
  geom_line(aes(time, hscore, color = factor(sigma2prior),group=interaction(sigma2prior,repl)))
#--------------------------------------------------------------------------------------------
# Check the partial log-evidence with training sample size m = 1 (RESCALED BY 1/t)
m = 1
partial_bayes = data.frame()
for (i in 1: length(sigma2prior_all)){
  for (r in 1:repl){
    logevidence = subset(results_all,sigma2prior== sigma2prior_all[i]&repl==r)$logevidence
    partiallogevidence = logevidence[(m+1):nobservations] - logevidence[1:m]
    partial_bayes = rbind(partial_bayes,data.frame(time = (m+1):nobservations, sigma2prior = sigma2prior_all[i], repl = r,
                                                   partiallogevidence = partiallogevidence))
  }
}
ggplot(partial_bayes) +
  geom_line(aes(time, -partiallogevidence, color = factor(sigma2prior),group=interaction(sigma2prior,repl)))

