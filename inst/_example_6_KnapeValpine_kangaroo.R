##################################################################################################
# Example 6: Kangaroo population models (Knape and Valpine, 2012)
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(wesanderson)
set.seed(19)

# Define data
timesteps = data_kangaroo["time",]
stepsize = 0.001
observations = data_kangaroo[1:2,]
nobservations = ncol(observations)

# Define hyperparameters
range_sigma = 10
range_tau = 10
range_r = 10
range_b = 10
# Define models
nb_models = 1
model = function(i){
  if (i==1){return(get_model_kangarooLogistic(timesteps,stepsize,range_sigma,range_tau,range_r,range_b))}
  if (i==2){return(get_model_kangarooExponential(timesteps,range_sigma,range_tau,range_r))}
  if (i==3){return(get_model_kangarooRandomwalk(timesteps,range_sigma,range_tau))}
}

# set algorithmic parameters
Nthetas = c(2^12,2^12,2^12)
Nxs = c(2^7,2^5,2^5)
algorithmic_parameters = list()
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------
##################################################################################################
##################################################################################################
##################################################################################################
#--------------------------------------------------------------------------------------------
### Compute logevidence and hscore for each model
results_all = data.frame()
post_all = list()
for (m in 1:nb_models){
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    algorithmic_parameters$Ntheta = Nthetas[m]
    algorithmic_parameters$Nx = Nxs[m]
    hscore(observations, model(m), algorithmic_parameters)
  }
  post = data.frame()
  for (r in 1:repl){
    results_all = rbind(results_all,data.frame(logevidence = results[[r]]$logevidence,
                                               hscore = results[[r]]$hscore,
                                               time = 1:nobservations,
                                               model = factor(m),
                                               repl = r))
    post = rbind(post,data.frame(t(results[[r]]$thetas_history[[nobservations+1]]),
                                 W = results[[r]]$normw_history[[nobservations+1]],
                                 repl = r))
  }
  post_all[[m]] = post
}
#--------------------------------------------------------------------------------------------
# Check posterior for each model
post_plot_all = list()
xlabels = list(xlab(expression(sigma)), xlab(expression(tau)), xlab("r"), xlab("b"))
colors = wes_palette("Darjeeling")[c(1,3,5)]
# Generate plots
for (m in 1:nb_models) {
  dimtheta = model(m)$dimtheta
  cnames = colnames(post_all[[m]])
  plot_post = list()
  for (k in 1:dimtheta) {
    local({k = k;
    cnames = cnames;
    plot_post[[k]] <<- ggplot(post_all[[m]]) +
      geom_density(aes_string(cnames[k],weight = "W", group = "repl"),col=colors[m],size=1,alpha=0.3) +
      theme(legend.position="none") + xlabels[[k]] + ylab("")})
  }
  post_plot_all[[m]] = plot_post
}
# display plots
do.call(grid.arrange,c(post_plot_all[[1]][c(3,4,1,2)], ncol = 2, nrow = 2))
do.call(grid.arrange,c(post_plot_all[[2]][c(3,1,2)], ncol = 3))
do.call(grid.arrange,c(post_plot_all[[3]], ncol = 2))

#--------------------------------------------------------------------------------------------
# Check the log-evidence
ggplot(results_all) +
  geom_line(aes(time, -logevidence/time, color = model,group=interaction(model,repl)), size = 1) +
  ylab("- log evidence / time") + guides(colour = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors)
#--------------------------------------------------------------------------------------------
# Check the h-score
ggplot(results_all) +
  geom_line(aes(time, hscore/time, color = model,group=interaction(model,repl)), size = 1) +
  ylab("Hyvarinen score / time") + guides(colour = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors)

