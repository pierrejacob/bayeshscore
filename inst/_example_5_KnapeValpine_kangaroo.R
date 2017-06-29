##################################################################################################
# Example 6: Kangaroo population models (Knap and Valpine, 2012)
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(wesanderson)
library(Hmisc)
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
nb_models = 3
model = function(i){
  if (i==1){return(get_model_kangarooLogistic(timesteps,stepsize,range_sigma,range_tau,range_r,range_b))}
  if (i==2){return(get_model_kangarooExponential(timesteps,range_sigma,range_tau,range_r))}
  if (i==3){return(get_model_kangarooRandomwalk(timesteps,range_sigma,range_tau))}
}
# set algorithmic parameters
Nthetas = c(2^14,2^12,2^12)
Nxs = c(2^5,2^5,2^5)

Nc = c(2^15,2^13,2^12)
Ncx = c(2^5,2^7,2^5)

algorithmic_parameters = list()
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$save = TRUE
algorithmic_parameters$store_last_byproducts = FALSE

algorithmic_parameters$nmoves = 2
algorithmic_parameters$reduce_variance = TRUE
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
#--------------------------------------------------------------------------------------------
# # keep track of progress in a log file
logfilename = "results.log"
writeLines(c(""), logfilename)
sink(file = logfilename, append = TRUE)
cat("Start time:",toString(Sys.time()),"\n")
#--------------------------------------------------------------------------------------------
for (m in 1:nb_models){
  cat("Model",toString(m)," starting time:",toString(Sys.time()),"\n")
  algorithmic_parameters$Ntheta = Nthetas[m]
  algorithmic_parameters$Nx = Nxs[m]
  algorithmic_parameters$Nc = Nc[m]
  algorithmic_parameters$Ncx = Ncx[m]
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    gc() # attempt at preserving RAM
    sink(file = logfilename, append = TRUE) # keep track of progress in a log file
    algorithmic_parameters$savefilename = paste("model_",toString(m),"_repl_",toString(i),".rds",sep="")
    hscore(observations, model(m), algorithmic_parameters)
  }
  post = data.frame()
  for (r in 1:repl){
    cat("Model",toString(m)," replication ",toString(r)," starting time:",toString(Sys.time()),"\n")
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
  savefilename = paste("final_result_kangaroo_model_",toString(m),"_",format(Sys.time(),"%Y-%m-%d_%H-%M-%S"),".rds",sep="")
  # saveRDS(list(results_all = results_all, post = post, model = m), file = savefilename)
  saveRDS(list(results = results, post = post, model = m), file = savefilename)
}
sink(file = NULL)
savefilename = paste("final_result_kangaroo","_",format(Sys.time(),"%Y-%m-%d_%H-%M-%S"),".rds",sep="")
saveRDS(list(results_all = results_all, post_all = post_all), file = savefilename)
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Generate plots for paper
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# plot posterior
post_plot_all = list()
colors = wes_palette("Darjeeling")[c(1,3,5)]
for (m in 1:nb_models){
  dimtheta = model(m)$dimtheta
  post = data.frame(thetas = unlist(post_all[[m]][,1:dimtheta]),
                    repl = rep(post_all[[m]]$repl,dimtheta),
                    w = rep(post_all[[m]]$W,dimtheta),
                    type = rep(factor(1:dimtheta),each=length(post_all[[m]]$repl)))
  if (m==1){levels(post$type) = c("r","b",expression(sigma),expression(tau)); nbcol = 2}
  if (m==2){levels(post$type) = c("r",expression(sigma),expression(tau)); nbcol = 3}
  if (m==3){levels(post$type) = c(expression(sigma),expression(tau)); nbcol = 2}

  local({m = m;
  post_plot_all[[m]] <<-ggplot(post) +
    geom_density(aes(thetas, weight = w, group = interaction(type,repl)),col=colors[m],size=1,alpha=0.3) +
    facet_wrap( ~ type, ncol=nbcol, scales="free",labeller=label_parsed) + xlab("") + ylab("Posterior density") +
    theme(strip.text.y = element_text(size = 12, colour = "black")) +
    theme(legend.text=element_text(size=12)) +
    theme(legend.title=element_text(size=12)) +
    theme(legend.position="none") +
    theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
    theme(axis.title.x=element_text(margin=margin(10,0,0,0)))})
}

post_plot_all[[1]]
# ggsave("example_5_kangaroos_posterior_model_1.png",plot = post_plot_all[[1]],dpi = 300,width = 10,height = 5)
post_plot_all[[2]]
# ggsave("example_5_kangaroos_posterior_model_2.png",plot = post_plot_all[[2]],dpi = 300,width = 10,height = 5)
post_plot_all[[3]]
# ggsave("example_5_kangaroos_posterior_model_3.png",plot = post_plot_all[[3]],dpi = 300,width = 10,height = 5)

#--------------------------------------------------------------------------------------------
colnames(results_all)[4] = "Model"
# plot log-evidence
ggplot(results_all) +
  geom_line(aes(time, -logevidence/time, color = Model,group=interaction(Model,repl)), size = 1) +
  ylab("- log evidence / time") + guides(colour = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave("example_5_kangaroos_logevidence.png",dpi = 300,width = 10,height = 5)


# plot Hyvarinen score
ggplot(results_all) +
  geom_line(aes(time, hscore/time, color = Model,group=interaction(Model,repl)), size = 1) +
  ylab("Prequential Hyvärinen score / time") + guides(colour = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors,"Model") +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave("example_5_kangaroos_preqhyvarinenscore.png",dpi = 300,width = 10,height = 5)

a= 0.95
ggplot(subset(results_all, time >= 0 & repl != 0),aes(time, -logevidence/time, color = Model)) +
  ylab("- log evidence / time") +
  # guides(shape = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  geom_line(aes(group=interaction(Model,repl)),linetype="dashed",alpha=0.5) +
  stat_summary(aes(group=Model,fill=Model),geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=a),alpha=0.3) +
  stat_summary(aes(group=Model,shape = Model),geom="point", fun.y=mean,size=2) +
  stat_summary(aes(group=Model),geom="line", fun.y=mean, size = 1) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave("example_5_kangaroos_logevidence_CI.png",dpi = 300,width = 10,height = 5)


ggplot(subset(results_all, time >= 0 & repl != 0),aes(time, hscore/time, color = Model)) +
  ylab("Prequential Hyvärinen score / time") +
  # guides(shape = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  geom_line(aes(group=interaction(Model,repl)),linetype="dashed",alpha=0.5) +
  stat_summary(aes(group=Model,fill=Model),geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=a),alpha=0.3) +
  stat_summary(aes(group=Model,shape = Model),geom="point", fun.y=mean,size=2) +
  stat_summary(aes(group=Model),geom="line", fun.y=mean, size = 1) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave("example_5_kangaroos_preqhyvarinenscore_CI.png",dpi = 300,width = 10,height = 5)

