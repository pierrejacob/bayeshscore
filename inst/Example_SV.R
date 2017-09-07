##################################################################################################
# Example: Levy-driven Stochastic Volatility models (Chopin, Jacob, Papaspiliopoulos, 2013)
# (see also Barndorff-Nielsen and Shephard, 2001)
##################################################################################################
# rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(wesanderson)
set.seed(19)
#--------------------------------------------------------------------------------------------
# Monitor progress in parallel via log file
setwd("C:/Users/shao/Desktop/HyvarinenSSM")
logfilename = "results.log"
writeLines(c(""), logfilename)
sink(logfilename, append = TRUE)
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$Nx = 2^7
algorithmic_parameters$Nx_max = 2^12
algorithmic_parameters$min_acceptance_rate = 0.10
algorithmic_parameters$nmoves = 1
algorithmic_parameters$verbose = TRUE
# use direct density-estimation
algorithmic_parameters$use_dde = TRUE
algorithmic_parameters$dde_options = list(Ny = 2^14, nb_steps = Inf)
# save intermediary results every "save_stepsize" observations
algorithmic_parameters$save = TRUE
algorithmic_parameters$save_stepsize = 100
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------
nobservations = 1000
timesteps = 1:nobservations
# define models
nb_models = 3
model = function(i){
  if(i==1){return(get_model_SVLevy_singlefactor(timesteps))}
  if(i==2){return(get_model_SVLevy_multifactor_noleverage(timesteps))}
  if(i==3){return(get_model_SVLevy_multifactor_withleverage(timesteps))}
}
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Case 1: simulated data from model 1 (single factor Levy-driven no leverage)
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# simulate observations
theta = c(0, 0, 0.5, 0.0625, 0.01)
observations = simulateData(model(1),theta,nobservations)$Y
# list of candidate models
models_to_run = c(1,2,3)
#--------------------------------------------------------------------------------------------
### Compute logevidence and hscore for each model
results_all = data.frame()
post_all = vector("list",length(models_to_run))
schedule = expand.grid(1:repl,models_to_run)
#--------------------------------------------------------------------------------------------
results = foreach(s=1:nrow(schedule),.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
  m = schedule[s,2]
  r = schedule[s,1]
  gc() # attempt to limit RAM usage
  sink(logfilename, append = TRUE) # Monitor progress in parallel via log file
  run = list()
  # for convenience, if error due to lack of RAM, re-run until completion (hacky but the while loop only kicks in when R runs out of RAM)
  while (length(run)==0){
    algorithmic_parameters$savefilename = paste("model_",toString(m),"_repl_",toString(r),"_",format(Sys.time(),"%Y-%m-%d_%H-%M-%S"),".rds",sep="")
    run =  tryCatch(hscore(observations, model(m), algorithmic_parameters),error = function(e){print(e);print("new attempt in progress ...");return (list())})
  }
  savefilename = paste("SV_run_model_",toString(m),"_repl_",toString(r),".rds",sep="")
  saveRDS(run, file = savefilename)
  run
}

for (m in models_to_run) {
  for (r in 1:repl) {
    s = (1:nrow(schedule))[(schedule[,1]==r) & (schedule[,2]==m)]
    # retrieve DDE estimates of the hscore
    hscoreDDE = results[[s]]$hscoreDDE
    if (algorithmic_parameters$dde_options$nb_steps < nobservations) {
      lastindex = algorithmic_parameters$dde_options$nb_steps
      hscoreDDE = c(hscoreDDE[1:lastindex], hscoreDDE[lastindex] + cumsum(diff(results[[s]]$hscore)[lastindex:(nobservations-1)]))
    }
    # format results into dataframes
    results_all = rbind(results_all,data.frame(logevidence = results[[s]]$logevidence,
                                               hscore = results[[s]]$hscore,
                                               hscoreDDE = hscoreDDE,
                                               time = 1:nobservations,
                                               model = factor(m),
                                               repl = r))
    post_all[[m]] = rbind(post_all[[m]],data.frame(t(results[[s]]$thetas),
                                                   W = results[[s]]$normw,
                                                   repl = r))

  }
  savefilename = paste("SV_model_",paste(1:m,sep="",collapse = "_"),".rds",sep="")
  saveRDS(list(post_all = post_all, results_all = results_all), file = savefilename)
}
#--------------------------------------------------------------------------------------------
# close log file
sink()
#--------------------------------------------------------------------------------------------
# Check posterior for each model
post_plot_all = list()
xlabels = list(xlab(expression(mu)), xlab(expression(beta)), xlab(expression(xi)),
               xlab(expression(omega^2)), xlab(expression(lambda[1])), xlab(expression(lambda[2])),
               xlab(expression(rho[1])), xlab(expression(rho[2])), xlab("w"))
colors = c(wes_palette("Darjeeling2")[c(2,3)],wes_palette("Darjeeling")[2])
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
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

g1 = do.call(grid.arrange,c(post_plot_all[[1]], ncol = 3, nrow = 3))
# ggsave("example_case3_SV_posterior_model_1.png",plot = g1,dpi = 300,width = 10, height = 5)

g2 = do.call(grid.arrange,c(post_plot_all[[2]], ncol = 3, nrow = 3))
# ggsave("example_case3_SV_posterior_model_2.png",plot = g2,dpi = 300,width = 10, height = 5)

g3 = do.call(grid.arrange,c(post_plot_all[[3]], ncol = 3, nrow = 3))
# ggsave("example_case3_SV_posterior_model_3.png",plot = g3,dpi = 300,width = 10, height = 5)


a = 0.95
# Check the log-evidence
ggplot(results_all,aes(time,  -logevidence, color = model)) +
  ylab("- log evidence") +
  # guides(shape = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  geom_line(aes(group=interaction(model,repl)),linetype="dashed",alpha=0.8) +
  # stat_summary(aes(group=model,fill=model),geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=a),alpha=0.3) +
  stat_summary(aes(group=model,shape = model),geom="point", fun.y=mean,size=2) +
  stat_summary(aes(group=model),geom="line", fun.y=mean, size = 1)
# ggsave("example_SV_logevidence.png",dpi = 300,width = 10, height = 5)

# Check the h-score DDE
ggplot(results_all,aes(time, hscoreDDE/time, color = model)) +
  ylab("Hyvarinen score / time") +
  # guides(shape = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  geom_line(aes(group=interaction(model,repl)),linetype="dashed",alpha=0.5) +
  # stat_summary(aes(group=model,fill=model),geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=a),alpha=0.3) +
  stat_summary(aes(group=model,shape = model),geom="point", fun.y=mean,size=2) +
  stat_summary(aes(group=model),geom="line", fun.y=mean, size = 1)
# ggsave("example_SV_hscore_rescaled.png",width = 10, height = 5,dpi = 300)

# Check the h-score DDE
ggplot(subset(results_all,time<50),aes(time, hscoreDDE, color = model)) +
  ylab("Hyvarinen score / time") +
  # guides(shape = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  geom_line(aes(group=interaction(model,repl)),linetype="dashed",alpha=0.5) +
  # stat_summary(aes(group=model,fill=model),geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=a),alpha=0.3) +
  stat_summary(aes(group=model,shape = model),geom="point", fun.y=mean,size=2) +
  stat_summary(aes(group=model),geom="line", fun.y=mean, size = 1)
# ggsave("example_SV_hscore.png",width = 10, height = 5,dpi = 300)
