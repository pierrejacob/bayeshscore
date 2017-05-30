##################################################################################################
# This implements example 3.2. in O'Hagan (1995).
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
set.seed(19)
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
#--------------------------------------------------------------------------------------------
# set hyperparameters
muprior = 0
sigma2prior = 100
kappa0 = 1
nu0 = 1
s02 = 1
# define models
nb_models = 3
model = function(i){
  if(i==1){return(get_model_iid_gaussian_unknown_mean(muprior,sigma2prior))} #iid N(theta1, 1)
  if(i==2){return(get_model_iid_gaussian_unknown_variance(nu0,s02))} #iid N(0, theta2)
  if(i==3){return(get_model_iid_gaussian(muprior,kappa0,nu0,s02))} #iid N(0, theta2)
}
#--------------------------------------------------------------------------------------------
nobservations = 100
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
##################################################################################################
# Case 1: DGP = N(0,5), Model 2 is correct
##################################################################################################
observations1 = matrix(rnorm(nobservations,0,sqrt(5)), nrow = 1)# observations (dimy by nobservations matrix)
#-----------------------------------
results1_all = data.frame()
post1_all = data.frame()
#-----------------------------------
# Compute logevidence and hscore
for (m in 1:nb_models){
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations1, model(m), algorithmic_parameters)
  }
  for (r in 1:repl){
    results1_all = rbind(results1_all,data.frame(logevidence = results[[r]]$logevidence,
                                                 hscore = results[[r]]$hscore,
                                                 time = 1:nobservations,
                                                 model = factor(m),
                                                 repl = factor(r)))
    post1_all = rbind(post1_all,data.frame(theta = c(results[[r]]$thetas_history[[nobservations+1]]),
                                           W = results[[r]]$normw_history[[nobservations+1]],
                                           model = factor(m),
                                           repl = factor(r)))
  }
}
##################################################################################################
# Case 2: DGP = N(1,1), Model 1 is correct
##################################################################################################
observations2 = matrix(rnorm(nobservations,1,1), nrow = 1)# observations (dimy by nobservations matrix)
#-----------------------------------
results2_all = data.frame()
post2_all = data.frame()
#-----------------------------------
# Compute logevidence and hscore
for (m in 1:nb_models){
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations2, model(m), algorithmic_parameters)
  }
  for (r in 1:repl){
    results2_all = rbind(results2_all,data.frame(logevidence = results[[r]]$logevidence,
                                                 hscore = results[[r]]$hscore,
                                                 time = 1:nobservations,
                                                 model = factor(m),
                                                 repl = factor(r)))
    post2_all = rbind(post2_all,data.frame(theta = c(results[[r]]$thetas_history[[nobservations+1]]),
                                           W = results[[r]]$normw_history[[nobservations+1]],
                                           model = factor(m),
                                           repl = factor(r)))
  }
}
##################################################################################################
# Case 3: DGP = N(2,3), both model 1 and 2 are misspecified
##################################################################################################
observations3 = matrix(rnorm(nobservations,2,sqrt(3)), nrow = 1)# observations (dimy by nobservations matrix)
#-----------------------------------
results3_all = data.frame()
post3_all = data.frame()
#-----------------------------------
# Compute logevidence and hscore
for (m in 1:nb_models){
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations3, model(m), algorithmic_parameters)
  }
  for (r in 1:repl){
    results3_all = rbind(results3_all,data.frame(logevidence = results[[r]]$logevidence,
                                                 hscore = results[[r]]$hscore,
                                                 time = 1:nobservations,
                                                 model = factor(m),
                                                 repl = factor(r)))
    post3_all = rbind(post3_all,data.frame(theta = c(results[[r]]$thetas_history[[nobservations+1]]),
                                           W = results[[r]]$normw_history[[nobservations+1]],
                                           model = factor(m),
                                           repl = factor(r)))
  }
}
##################################################################################################
# Case 4: DGP = N(0,1), both model 1 and 2 are correct
##################################################################################################
observations4 = matrix(rnorm(nobservations,0,1), nrow = 1)# observations (dimy by nobservations matrix)
#-----------------------------------
results4_all = data.frame()
post4_all = data.frame()
#-----------------------------------
# Compute logevidence and hscore
for (m in 1:nb_models){
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations4, model(m), algorithmic_parameters)
  }
  for (r in 1:repl){
    results4_all = rbind(results4_all,data.frame(logevidence = results[[r]]$logevidence,
                                                 hscore = results[[r]]$hscore,
                                                 time = 1:nobservations,
                                                 model = factor(m),
                                                 repl = factor(r)))
    post4_all = rbind(post4_all,data.frame(theta = c(results[[r]]$thetas_history[[nobservations+1]]),
                                           W = results[[r]]$normw_history[[nobservations+1]],
                                           model = factor(m),
                                           repl = factor(r)))
  }
}
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
#--------------------------------------------------------------------------------------------
# Checking h score
#--------------------------------------------------------------------------------------------
results_all = list(results1_all, results2_all, results3_all, results4_all)
plot_hscore = list()
for (i in 1:4){
  local({
    i = i;
    results_all = results_all;
    plot_hscore[[i]] <<-ggplot() +
      geom_line(data=results_all[[i]],aes(time,hscore/time,group=interaction(repl,model),col=model))
  })
}
do.call(grid.arrange,c(plot_hscore, ncol = 2))

