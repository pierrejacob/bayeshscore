##################################################################################################
# Sanity check ARMA(1,0) vs AR(1)
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
  if (i==1){return(get_model_ARMA(1,0,nu0, sigma02))} # ARMA(1,0)
  if (i==2){return(get_model_AR1(nu0, sigma02))} #AR(1)
}

# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R

# generate some observations
nobservations = 20
true_theta = c(0.5,2)
observations = simulateData(model(1),true_theta,nobservations)$Y
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------
results_all = data.frame()
post_all = data.frame()
#--------------------------------------------------------------------------------------------
### Compute logevidence and hscore for increasing vagueness
for (m in 1:nb_models){
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations, model(m), algorithmic_parameters)
  }
  for (r in 1:repl){
    results_all = rbind(results_all,data.frame(logevidence = results[[r]]$logevidence,
                                               hscore = results[[r]]$hscore,
                                               time = 1:nobservations,
                                               model = factor(m),
                                               repl = r))
    post_all = rbind(post_all,data.frame(phi1 = c(results[[r]]$thetas_history[[nobservations+1]][1,]),
                                         sigma2 = c(results[[r]]$thetas_history[[nobservations+1]][2,]),
                                         W = results[[r]]$normw_history[[nobservations+1]],
                                         model = factor(m),
                                         repl = factor(r)))
  }
}
#--------------------------------------------------------------------------------------------
### Check posterior
g1 = ggplot(post_all) +
  geom_density(aes(phi1, weight = W, color = model, group = interaction(repl,model)),alpha=0.3,size=1) +
  theme(legend.position="none") + xlab(expression(varphi[1])) + ylab("")
g2 = ggplot(post_all) +
  geom_density(aes(sigma2, weight = W, color = model, group = interaction(repl,model)),alpha=0.3,size=1) +
  xlab(expression(sigma^2)) + ylab("")
grid.arrange(g1,g2,widths=c(1,1.15))
#--------------------------------------------------------------------------------------------
### Check logevidence
ggplot(results_all) +
  geom_line(aes(time,-logevidence/time,color=model,group=interaction(model,repl)))
#--------------------------------------------------------------------------------------------
### Check H-score
ggplot(results_all) +
  geom_line(aes(time,hscore/time,color=model,group=interaction(model,repl)))
