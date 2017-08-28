##################################################################################################
# Example 4: AR(1) vs MA(1)
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
all_models = list()
all_models[[1]] = list(p = 1, q = 0)
all_models[[2]] = list(p = 0, q = 1)
model = function(i){
  return(get_model_ARMA(all_models[[i]]$p,all_models[[i]]$q,nu0, sigma02))
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
algorithmic_parameters$Nx = 2^7
algorithmic_parameters$verbose = TRUE

algorithmic_parameters$nmoves = 5
# algorithmic_parameters$reduce_variance = TRUE
# algorithmic_parameters$Nc = 2^12
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------
nobservations = 1000

##################################################################################################
# Case 1: true model = AR(1)
##################################################################################################
true_model = 1
# true_theta = get_stationaryparameters(all_models[[true_model]]$p,all_models[[true_model]]$q)
true_theta = c(0.5,1)
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
# true_theta = get_stationaryparameters(all_models[[true_model]]$p,all_models[[true_model]]$q)
true_theta = c(0.5,1)
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
    logbayesfactors = rbind(logbayesfactors,data.frame(case = factor(i),
                                                       time = 1:nobservations,
                                                       repl = r,
                                                       logbayesfactor = logbayes_factor,
                                                       type = factor(paste("Case ",toString(i)))))
    h_factor = subset(results,model==2&repl==r)$hscore - subset(results,model==1&repl==r)$hscore
    h_factors = rbind(h_factors,data.frame(case = factor(i),
                                           time = 1:nobservations,
                                           repl = r,
                                           hfactor = h_factor,
                                           type = factor(paste("Case ",toString(i)))))
    # local({i = i;
    # BF_plots[[i]] <<- ggplot(subset(logbayesfactors, case==i)) +
    #   geom_line(aes(time, logbayesfactor, color = case, group = repl)) +
    #   geom_hline(yintercept = 0,linetype="dotted",size=1) +
    #   ylab("log Bayes factor");
    # HF_plots[[i]] <<- ggplot(subset(h_factors, case==i)) +
    #   geom_line(aes(time, hfactor, color = case, group = repl)) +
    #   geom_hline(yintercept = 0,linetype="dotted",size=1) +
    #   ylab("H factor")
    # })
  }
}
# # Plot log Bayes factor
# # left, right = case 1, 2
# do.call(grid.arrange,c(BF_plots, ncol = 2))
# # Plot H factor
# # left, right = case 1, 2
# do.call(grid.arrange,c(HF_plots, ncol = 2))
##################################################################################################
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Generate plots for paper
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Bayes factor
ggplot(logbayesfactors) +
  geom_line(aes(time, logbayesfactor, color = case, group = interaction(case,repl))) +
  geom_hline(yintercept = 0,linetype="dotted",size=1) +
  ylab("log Bayes factor  [1 vs 2]") + facet_grid(. ~ type) + xlab("Number of observations") +
  # guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))

# ggsave("example_4_AR1_MA1_log_BF_1_vs_2.png",dpi = 300,width = 10,height = 5)

case_label <- list(
  'Case 1'=expression(paste("Case 1: ",M[1]," is well-specified",sep="")),
  'Case 2'=expression(paste("Case 2: ",M[2]," is well-specified",sep=""))
)
case_labeller <- function(variable,value){
  return(case_label[value])
}
# Hyvarinen factor
ggplot(h_factors) +
  geom_line(aes(time, hfactor, color = case, group = interaction(case,repl))) +
  geom_hline(yintercept = 0,linetype="dotted",size=1) +
  xlab("Number of observations") +
  ylab("Hyvrärinen factor  [1 vs 2]") +
  facet_grid(. ~ type, labeller = case_labeller) +
  # guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))

# ggsave("example_4_AR1_MA1_Hyvarinen_factor_1_vs_2.png",dpi = 300,width = 10,height = 5)
