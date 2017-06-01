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

# this function finds parameters (used to generate data) that guarantee stationarity
get_stationaryparameters = function(p,q){
  AR_coeffs = runif(p,-1,1)
  MA_coeffs = runif(q,-1,1)
  AR_roots = polyroot(c(1,-AR_coeffs))
  MA_roots = polyroot(c(1,-MA_coeffs))
  accept = ((sum(abs(AR_roots)<=1))==0)&&((sum(abs(MA_roots)<=1))==0)&&(length(intersect(AR_roots,MA_roots))==0)
  while (!accept){
    AR_coeffs = runif(p,-1,1)
    MA_coeffs = runif(q,-1,1)
    accept = ((sum(abs(polyroot(c(1,-AR_coeffs)))<=1))==0)&&((sum(abs(polyroot(c(1,-MA_coeffs)))<=1))==0)
  }
  return (c(AR_coeffs,MA_coeffs,1))
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
nobservations = 50

##################################################################################################
# Case 1: true model = AR(1)
##################################################################################################
true_model = 1
true_theta = get_stationaryparameters(1,0)
observations1 = simulateData(model(true_model),true_theta,nobservations)$Y
# observations in a matrix of dimensions dimy x nobservations
#--------------------------------------------------------------------------------------------
results_all1 = data.frame()
### Compute logevidence and hscore for each model
for (m in 1:nb_models){
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations1, model(m), algorithmic_parameters)
  }
  for (r in 1:repl){
    results_all1 = rbind(results_all1,data.frame(logevidence = results[[r]]$logevidence,
                                                 hscore = results[[r]]$hscore,
                                                 time = 1:nobservations,
                                                 model = factor(m),
                                                 repl = r))
  }
}
### Plot log evidence
ggplot(results_all1) +
  geom_line(aes(time,-logevidence/time,color=model,group=interaction(model,repl)))
### Plot H score
ggplot(results_all1) +
  geom_line(aes(time,hscore/time,color=model,group=interaction(model,repl)))

##################################################################################################
# Case 2: true model = MA(1)
##################################################################################################
true_model = 2
true_theta = get_stationaryparameters(0,1)
observations2 = simulateData(model(true_model),true_theta,nobservations)$Y
# observations in a matrix of dimensions dimy x nobservations
#--------------------------------------------------------------------------------------------
results_all2 = data.frame()
### Compute logevidence and hscore for each model
for (m in 1:nb_models){
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations2, model(m), algorithmic_parameters)
  }
  for (r in 1:repl){
    results_all2 = rbind(results_all2,data.frame(logevidence = results[[r]]$logevidence,
                                                 hscore = results[[r]]$hscore,
                                                 time = 1:nobservations,
                                                 model = factor(m),
                                                 repl = r))
  }
}
### Plot log evidence
ggplot(results_all2) +
  geom_line(aes(time,-logevidence/time,color=model,group=interaction(model,repl)))
### Plot H score
ggplot(results_all2) +
  geom_line(aes(time,hscore/time,color=model,group=interaction(model,repl)))


##################################################################################################
# Compare model 1 and model 2 in each case
##################################################################################################
results_all = list(results_all1,results_all2)
logbayesfactors = data.frame()
h_factors = data.frame()
BF_plots = list()
HF_plots = list()
for (r in 1:repl) {
  for (i in 1:nb_models) {
    results = results_all[[i]]
    logbayes_factor = subset(results,model==1&repl==r)$logevidence - subset(results,model==2&repl==r)$logevidence
    logbayesfactors = rbind(logbayesfactors,data.frame(case = factor(i), time = 1:nobservations, repl = r, logbayesfactor = logbayes_factor))
    h_factor = subset(results,model==2&repl==r)$hscore - subset(results,model==1&repl==r)$hscore
    h_factors = rbind(h_factors,data.frame(case = factor(i), time = 1:nobservations, repl = r, hfactor = h_factor))
    local({i = i;
    BF_plots[[i]] <<- ggplot(subset(logbayesfactors, case==i)) +
      geom_line(aes(time, logbayesfactor, color = case, group = repl)) +
      geom_hline(yintercept = 0,linetype="dotted",size=1) +
      ylab("log Bayes factor");
    HF_plots[[i]] <<- ggplot(subset(h_factors, case==i)) +
      geom_line(aes(time, hfactor, color = case, group = repl)) +
      geom_hline(yintercept = 0,linetype="dotted",size=1) +
      ylab("H factor")
    })
  }
}
# Plot log Bayes factor
# left, right = case 1, 2
do.call(grid.arrange,c(BF_plots, ncol = 2))
# Plot H factor
# left, right = case 1, 2
do.call(grid.arrange,c(HF_plots, ncol = 2))

