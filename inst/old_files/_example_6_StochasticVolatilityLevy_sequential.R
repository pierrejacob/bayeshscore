##################################################################################################
# Example 7: Levy-driven Stochastic Volatility models (Chopin, Jacob, Papaspiliopoulos, 2013)
# (see also Barndorff-Nielsen and Shephard, 2001)
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(wesanderson)
set.seed(19)
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$Nx = 2^7
algorithmic_parameters$Nxmax = 2^12
algorithmic_parameters$nmoves = 5
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$save = FALSE

algorithmic_parameters$reduce_variance = TRUE
algorithmic_parameters$Nc = 2^12
algorithmic_parameters$Ncx = 2^10

algorithmic_parameters$use_kde = TRUE
algorithmic_parameters$kde_opt = list(Ny = 10^4, nb_steps = nobservations)
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Case 1: simulated data from model 1 (single factor Levy-driven)
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# simulate observations
nobservations = 100
timesteps = 1:nobservations
theta = c(0, 0, 0.5, 0.0625, 0.01)
observations = simulateData(get_model_SVLevy_singlefactor(timesteps),theta,nobservations)$Y
#--------------------------------------------------------------------------------------------
# define models
nb_models = 3
model = function(i){
  if(i==1){return(get_model_SVLevy_singlefactor(timesteps))}
  if(i==2){return(get_model_SVLevy_multifactor_noleverage(timesteps))}
  if(i==3){return(get_model_SVLevy_multifactor_withleverage(timesteps))}
}
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
models_to_run = c(2,3)
#--------------------------------------------------------------------------------------------
### Compute logevidence and hscore for each model
results_all = data.frame()
post_all = list()
#--------------------------------------------------------------------------------------------
for (m in models_to_run){
  gc() # attempt to limit RAM usage
  cat("Model ",toString(m)," started at:", toString(Sys.time()),"\n")
  post = data.frame()
  for (r in 1:repl){
    results = hscore(observations, model(m), algorithmic_parameters)
    results_all = rbind(results_all,data.frame(logevidence = results[[r]]$logevidence,
                                               hscore = results[[r]]$hscore,
                                               time = 1:nobservations,
                                               model = factor(m),
                                               repl = r))
    post = rbind(post,data.frame(t(results[[r]]$thetas),
                                 W = results[[r]]$normw,
                                 repl = r))
  }
  post_all[[m]] = post
}
#--------------------------------------------------------------------------------------------
# Check posterior for each model
post_plot_all = list()
xlabels = list(xlab(expression(mu)), xlab(expression(beta)), xlab(expression(xi)),
               xlab(expression(omega^2)), xlab(expression(lambda[1])), xlab(expression(lambda[2])),
               xlab(expression(rho[1])), xlab(expression(rho[2])), xlab("w"))
colors = wes_palette("Darjeeling")[c(1,3,5)]
# Generate plots
for (m in models_to_run) {
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
for (m in models_to_run) {
  do.call(grid.arrange,c(post_plot_all[[m]], ncol = 3, nrow = 3))
}

#--------------------------------------------------------------------------------------------
# Check the log-evidence
ggplot(results_all) +
  geom_line(aes(time, -logevidence/time, color = model,group=interaction(model,repl)), size = 1) +
  ylab("- log evidence / time") + guides(colour = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors)
#--------------------------------------------------------------------------------------------
# Check the h-score
ggplot(subset(results_all, time >= 1)) +
  geom_line(aes(time, hscore/time, color = model,group=interaction(model,repl)), size = 1) +
  ylab("Hyvarinen score / time") + guides(colour = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors)
# #--------------------------------------------------------------------------------------------
# # Checking h factor
# #--------------------------------------------------------------------------------------------
# h_factors = data.frame()
# plot_hfactor = list()
# for (m in 1:(nb_models-1)){
#   for (r in 1:repl) {
#     case = factor(paste(toString(m+1),"/ 1"))
#     h_factor = subset(results_all,model==(m+1)&repl==r)$hscore-subset(results_all,model==1&repl==r)$hscore
#     h_factors = rbind(h_factors,data.frame(time = 1:nobservations,
#                                            repl = r,
#                                            hfactor = h_factor,
#                                            model = case))
#   }
#   local({m = m;
#   case = case;
#   plot_hfactor[[m]] <<- ggplot(subset(h_factors, case = case)) +
#     geom_line(aes(time, hfactor, group = repl)) +
#     ylab("H factor")
#   })
# }
# # Checking H-factor
# # top-left, top-right, bottom-left, bottom-right = case 1, 2, 3, 4
# do.call(grid.arrange,c(plot_hfactor, ncol = 2))
