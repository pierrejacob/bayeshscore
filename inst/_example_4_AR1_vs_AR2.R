##################################################################################################
# Example 4: AR(1) vs AR(2)
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
  if (i==1){return(get_model_AR1(nu0, sigma02))}
  if (i==2){return(get_model_AR2(nu0, sigma02))}
}

# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^12
algorithmic_parameters$verbose = TRUE
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------


##################################################################################################
# Case 1: true model = AR(1)
##################################################################################################
nobservations = 100
true_model = 1
true_theta = c(0.5,2)
observations1 = simulateData(model(true_model),true_theta,nobservations)
# observations in a matrix of dimensions dimy x nobservations
#--------------------------------------------------------------------------------------------
results_all1 = data.frame()
#--------------------------------------------------------------------------------------------
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
#--------------------------------------------------------------------------------------------
# Check the log-evidence
ggplot(results_all1) +
  geom_line(aes(time, -logevidence/time, color = model,group=interaction(model,repl))) +
  ylab("- log evidence") + guides(colour = guide_legend(override.aes = list(size=2)))
#--------------------------------------------------------------------------------------------
# Check the h-score
ggplot(results_all1) +
  geom_line(aes(time, hscore/time, color = model,group=interaction(model,repl))) +
  ylab("Hyvarinen score") + guides(colour = guide_legend(override.aes = list(size=2)))


##################################################################################################
# Case 2: true model = AR(2)
##################################################################################################
true_model = 2
true_theta = c(0.5,0.25,2)
observations2 = simulateData(model(true_model),true_theta,nobservations)
# observations in a matrix of dimensions dimy x nobservations
#--------------------------------------------------------------------------------------------
results_all2 = data.frame()
#--------------------------------------------------------------------------------------------
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
#--------------------------------------------------------------------------------------------
# Check the log-evidence
ggplot(results_all2) +
  geom_line(aes(time, -logevidence/time, color = model,group=interaction(model,repl))) +
  ylab("- log evidence") + guides(colour = guide_legend(override.aes = list(size=2)))
#--------------------------------------------------------------------------------------------
# Check the h-score
ggplot(results_all2) +
  geom_line(aes(time, hscore/time, color = model,group=interaction(model,repl))) +
  ylab("Hyvarinen score") + guides(colour = guide_legend(override.aes = list(size=2)))


##################################################################################################
# Plot results
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
