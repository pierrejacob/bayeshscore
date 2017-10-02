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
registerDoParallel(cores=5) #number of workers in parallel
models_to_run = c(1,2,3)
#--------------------------------------------------------------------------------------------
### Compute logevidence and hscore for each model
results_all = data.frame()
post_all = list()
#--------------------------------------------------------------------------------------------
# Monitor progress in parallel via log file
logfilename = "results.log"
writeLines(c(""), logfilename)
sink(logfilename, append = TRUE)
#--------------------------------------------------------------------------------------------
for (m in models_to_run){
  gc() # attempt to limit RAM usage
  cat("Model ",toString(m)," started at:", toString(Sys.time()))
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    gc() # attempt to limit RAM usage
    sink(logfilename, append = TRUE) # Monitor progress in parallel via log file
    # algorithmic_parameters$savefilename = paste("model_",toString(m),"_repl_",toString(i),"_",format(Sys.time(),"%Y-%m-%d_%H-%M-%S"),".rds",sep="")
    hscore(observations, model(m), algorithmic_parameters)
  }
  post = data.frame()
  for (r in 1:repl){
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
# for (m in 1:nb_models){
#   cat("Model ",toString(m)," started at:", toString(Sys.time()))
#   results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
#     sink(logfilename, append = TRUE) # Monitor progress in parallel via log file
#     algorithmic_parameters$savefilename = paste("model_",toString(m),"_repl_",toString(i),".rds",sep="")
#     hscore(observations, model(m), algorithmic_parameters)
#   }
#   post = data.frame()
#   for (r in 1:repl){
#     results_all = rbind(results_all,data.frame(logevidence = results[[r]]$logevidence,
#                                                hscore = results[[r]]$hscore,
#                                                time = 1:nobservations,
#                                                model = factor(m),
#                                                repl = r))
#     post = rbind(post,data.frame(t(results[[r]]$thetas),
#                                  W = results[[r]]$normw,
#                                  repl = r))
#   }
#   post_all[[m]] = post
# }
#--------------------------------------------------------------------------------------------
# close log file
sink()


#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Retireve results from previous runs
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
### Compute logevidence and hscore for each model
setwd("C:/Users/jianfushao/Documents/Harvard/MyHarvard/_Research_/Model Selection for SSM (Pierre and Jie)/Simulation/Example 6 - Stochastic Volatility (Chopin et al., 2013)")
results_all = data.frame()
post_all = list()
repl = 10 #number of replications
models_to_run = c(1,2)
#--------------------------------------------------------------------------------------------
for (m in models_to_run){
  post = data.frame()
  for (r in 1:repl){
    results = readRDS(paste("SV_model_",toString(m),"_repl_",toString(r),".rds",sep=""))
    results_all = rbind(results_all,data.frame(logevidence = results$logevidence,
                                               hscore = results$hscore,
                                               hscoreKDE = results$hscoreKDE,
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
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Generate plots
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Check posterior for each model
post_plot_all = list()
xlabels = list(xlab(expression(mu)), xlab(expression(beta)), xlab(expression(xi)),
               xlab(expression(omega^2)), xlab(expression(lambda[1])), xlab(expression(lambda[2])),
               xlab(expression(rho[1])), xlab(expression(rho[2])), xlab("w"))
colors = c(wes_palette("Darjeeling2")[2],wes_palette("FantasticFox")[c(4,5)])

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
results_to_plot = subset(results_all, time>=0 & time<=50 & repl!=0)
# Check the log-evidence
ggplot(results_to_plot) +
  geom_line(aes(time, -logevidence/time, color = model,group=interaction(model,repl)), size = 1) +
  ylab("- log evidence / time") + guides(colour = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors)
# ggsave("example_6_SV_logevidence.png",dpi = 300,width = 10,height = 5)
#--------------------------------------------------------------------------------------------
# # Check the h-score
# ggplot(subset(results_all, time >= 1 & repl != 0)) +
#   geom_line(aes(time, hscore/time, color = model,group=interaction(model,repl)), size = 1) +
#   ylab("Hyvarinen score / time") + guides(colour = guide_legend(override.aes = list(size=2))) +
#   scale_color_manual(values = colors)
#--------------------------------------------------------------------------------------------
# Check the h-score KDE
ggplot(results_to_plot) +
  geom_line(aes(time, hscoreKDE/time, color = model,group=interaction(model,repl)), size = 1) +
  ylab("Prequenial Hyvärinen score / time") + guides(colour = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors)
# ggsave("example_6_SV_preqhyvarinenscore.png",dpi = 300,width = 10,height = 5)


# Check the log-evidence
a = 0.95

ggplot(results_to_plot,aes(time,-logevidence/time,color=model)) +
  ylab("- log evidence / time") +
  # guides(shape = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  geom_line(aes(group=interaction(model,repl)),linetype="dashed",alpha=0.5) +
  stat_summary(aes(group=model,fill=model),geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=a),alpha=0.3) +
  stat_summary(aes(group=model,shape = model),geom="point", fun.y=mean,size=2) +
  stat_summary(aes(group=model),geom="line", fun.y=mean, size = 1)
# ggsave("example_6_SV_logevidence_CI.png",dpi = 300,width = 10,height = 5)

# Check the h-score KDE
ggplot(results_to_plot,aes(time,hscoreKDE/time,color=model)) +
  ylab("Prequenial Hyvärinen score / time") +
  # guides(shape = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  geom_line(aes(group=interaction(model,repl)),linetype="dashed",alpha=0.5) +
  stat_summary(aes(group=model,fill=model),geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=a),alpha=0.3) +
  stat_summary(aes(group=model,shape = model),geom="point", fun.y=mean,size=2) +
  stat_summary(aes(group=model),geom="line", fun.y=mean, size = 1)
# ggsave("example_6_SV_preqhyvarinenscore_CI.png",dpi = 300,width = 10,height = 5)


colors = wes_palette("Rushmore")[5]
hfactor1v2 = data.frame(time=subset(results_all,model==1)$time,
                        repl = subset(results_all,model==1)$repl,
                        hfactor = subset(results_all,model==1)$hscoreKDE-subset(results_all,model==2)$hscoreKDE)
ggplot(hfactor1v2,aes(time,hfactor)) +
  ylab("Hyvärinen factor [1 vs 2]") +
  # guides(shape = guide_legend(override.aes = list(size=2))) +
  # scale_color_manual(values = colors) +
  # scale_fill_manual(values = colors) +
  stat_summary(geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=a),alpha=0.3,fill=colors) +
  stat_summary(geom="point", fun.y=mean,size=2,color = colors) +
  stat_summary(geom="line", fun.y=mean, size = 1,color = colors) +
  geom_line(aes(group=repl),linetype="dashed",alpha=0.5,color = colors) +
  geom_hline(yintercept = 0, alpha = 0.3)

