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
algorithmic_parameters$Ntheta = 2^7
algorithmic_parameters$Nx = 2^7
algorithmic_parameters$nmoves = 5
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$save = FALSE
algorithmic_parameters$store_X_history = FALSE
algorithmic_parameters$store_thetas_history = FALSE
algorithmic_parameters$store_last_X = FALSE
algorithmic_parameters$store_last_thetas = TRUE
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Case 1: simulated data from model 1 (single factor Levy-driven)
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# simulate observations
nobservations = 5
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
#--------------------------------------------------------------------------------------------
### Compute logevidence and hscore for each model
results_all = data.frame()
post_all = list()
#--------------------------------------------------------------------------------------------
for (m in 1:nb_models){
  post = data.frame()
  for (r in 1:repl){
    cat("Model ",toString(m)," repl ",toString(r)," started at:", toString(Sys.time()),"\n")
    results = hscore(observations, model(m), algorithmic_parameters)
    results_all = rbind(results_all,data.frame(logevidence = results$logevidence,
                                               hscore = results$hscore,
                                               time = 1:nobservations,
                                               model = factor(m),
                                               repl = r))
    post = rbind(post,data.frame(t(results$thetas),
                                 W = results$normw,
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
do.call(grid.arrange,c(post_plot_all[[1]], ncol = 3, nrow = 3))
do.call(grid.arrange,c(post_plot_all[[2]], ncol = 3, nrow = 3))
do.call(grid.arrange,c(post_plot_all[[3]], ncol = 3, nrow = 3))
#--------------------------------------------------------------------------------------------
# Check the log-evidence
ggplot(results_all) +
  geom_line(aes(time, -logevidence/time, color = model,group=interaction(model,repl)), size = 1) +
  ylab("- log evidence / time") + guides(colour = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors)
#--------------------------------------------------------------------------------------------
# Check the h-score
ggplot(subset(results_all)) +
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
