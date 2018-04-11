##################################################################################################
# Example - Kangaroo population models (Knape and Valpine, 2012)
##################################################################################################
library(bayeshscore)
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

# Plot kangaroo data
kangaroo.df = data.frame(t(data_kangaroo))
ggplot(kangaroo.df, aes(x = time)) +
  geom_segment(aes(x = time, y = y1, xend = time, yend = y2),linetype="dashed",colour="black") +
  geom_point(aes(y=y1),size=2,colour = wes_palette("Darjeeling2")[2]) +
  geom_point(aes(y=y2),size=2,colour = wes_palette("Darjeeling")[1]) +
  ylab("Number of red kangaroos\n") +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0))) +
  xlab("\n Time (years)")
# ggsave("example_5_kangaroos_data.png",dpi = 300,width = 10,height = 5)

# Define hyperparameters
range_sigma = 10
range_tau = 10
range_r = 10
range_b = 10
# Define models
model = function(i){
  if (i==1){return(get_model_kangarooLogistic(timesteps,stepsize,range_sigma,range_tau,range_r,range_b))}
  if (i==2){return(get_model_kangarooExponential(timesteps,range_sigma,range_tau,range_r))}
  if (i==3){return(get_model_kangarooRandomwalk(timesteps,range_sigma,range_tau))}
}
models_to_run = c(1,2,3)
# set algorithmic parameters
Nthetas = c(2^14,2^14,2^14)
Nxs = c(2^5,2^5,2^5)

algorithmic_parameters = list()
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$save = TRUE
algorithmic_parameters$save_stepsize = 5
algorithmic_parameters$store_last_byproducts = FALSE
algorithmic_parameters$nmoves = 2
algorithmic_parameters$discrete_diff_type = "central"
algorithmic_parameters$resampling = function(normw) ssp_resampling_n(normw, runif(length(normw)))

# algorithmic_parameters$reduce_variance = TRUE
# Nc = c(2^15,2^13,2^12)
# Ncx = c(2^5,2^7,2^5)
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
for (m in models_to_run){
  cat("Model",toString(m)," starting time:",toString(Sys.time()),"\n")
  algorithmic_parameters$Ntheta = Nthetas[m]
  algorithmic_parameters$Nx = Nxs[m]
  algorithmic_parameters$Nc = Nc[m]
  algorithmic_parameters$Ncx = Ncx[m]
  results = foreach(i=1:repl,.packages=c('bayeshscore'),.verbose = TRUE) %dorng% {
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
                                               Model = factor(m),
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
axis_titlesize = 24
axis_ticktextsize = 20
colors = wes_palette("Darjeeling")[c(1,3,5)]
for (m in models_to_run){
  dimtheta = model(m)$dimtheta
  post = data.frame(thetas = unlist(post_all[[m]][,1:dimtheta]),
                    repl = rep(post_all[[m]]$repl,dimtheta),
                    w = rep(post_all[[m]]$W,dimtheta),
                    type = rep(factor(1:dimtheta),each=length(post_all[[m]]$repl)))
  if (m==1){levels(post$type) = c(expression(sigma),expression(tau),"r","b"); nbcol = 4}
  if (m==2){levels(post$type) = c(expression(sigma),expression(tau),"r"); nbcol = 3}
  if (m==3){levels(post$type) = c(expression(sigma),expression(tau)); nbcol = 2}

  local({m = m;
  post_plot_all[[m]] <<-ggplot(post) +
    geom_density(aes(thetas, weight = w, group = interaction(type,repl)),adjust = 1.25,col=colors[m],size=1,alpha=0.8) +
    facet_wrap( ~ type, ncol=nbcol, scales="free",labeller=label_parsed) + xlab("") + ylab("") +
    theme(axis.text.x = element_text(size = axis_ticktextsize),
          axis.text.y = element_text(size = axis_ticktextsize),
          axis.title.x = element_text(size = axis_titlesize, margin=margin(20,0,0,0)),
          axis.title.y = element_text(size = axis_titlesize, angle = 90, margin = margin(0,20,0,0)),
          strip.text.y = element_text(size = axis_titlesize, colour = "black"),
          strip.text.x = element_text(size = axis_titlesize, colour = "black"),
          strip.background = element_rect(fill="gray88"),
          panel.background = element_rect(fill="gray95",linetype = "solid", colour="white"),
          legend.position = "none")})
}

post_plot_all[[1]]
post_plot_all[[2]]
post_plot_all[[3]]
# ggsave("example_kangaroos_posterior_model_1_24_by_6.png",plot = post_plot_all[[1]],dpi = 300,width = 24,height = 6)
# ggsave("example_kangaroos_posterior_model_2_24_by_6.png",plot = post_plot_all[[2]],dpi = 300,width = 24,height = 6)
# ggsave("example_kangaroos_posterior_model_3_24_by_6.png",plot = post_plot_all[[3]],dpi = 300,width = 24,height = 6)
# ggsave("example_kangaroos_posterior_model_1_24_by_6.pdf",plot = post_plot_all[[1]],dpi = 300,width = 24,height = 6)
# ggsave("example_kangaroos_posterior_model_2_24_by_6.pdf",plot = post_plot_all[[2]],dpi = 300,width = 24,height = 6)
# ggsave("example_kangaroos_posterior_model_3_24_by_6.pdf",plot = post_plot_all[[3]],dpi = 300,width = 24,height = 6)

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
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
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
labels.df = data.frame(x = rep(44,3), y = c(-0.0040,-0.00445,-0.0048),
                       text = c("Model 1","Model 2","Model 3"),
                       type = factor(c("Model 1","Model 2","Model 3")))
ggplot() +
  ylab("Prequential Hyvärinen score") +
  xlab("Time") +
  xlim(0,45) +
  # geom_label(aes(rep(44,3),c(-0.0042,-0.0047,-0.0052),label = c("Model 1","Model 2","Model 3"), fontface = "bold"),
  # fill=colors, color = "black",alpha=0.5) +
  geom_label(data = labels.df, aes(x,y,label = text,color=type), color = colors, fontface = "bold") +
  theme(legend.position="none") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  # geom_line(data = results_all,aes(time, hscore, color = Model,group=interaction(Model,repl),linetype=Model)) +
  # scale_linetype_manual(values = c("dotted","solid","longdash")) +
  # stat_summary(data = results_all,aes(time, hscore, color = Model,group=Model,fill=Model),geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=0.95),alpha=0.3) +
  geom_line(data = results_all,aes(time, hscore, color = Model, group=interaction(Model,repl)),alpha=0.7) +
  stat_summary(data = results_all,aes(time, hscore, color = Model,group=Model,shape = Model),geom="point", fun.y=mean,size=2.5) +
  stat_summary(data = results_all,aes(time, hscore, color = Model,group=Model),geom="line", fun.y=mean, size = 1) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave("example_kangaroos_preqhyvarinenscore.png",dpi = 300,width = 10,height = 5)
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#### Artificially put all the things to be plotted into a single dataframe in order to use facets
results.df = results_all
results.df[1] = -results.df[1] # Turn the log-evidence into NEGATIVE log-evidence
colnames(results.df)[1:2] = rep("score")
results.df = rbind(results.df[,-2],results.df[,-1])
results.df = cbind(results.df,rep(factor(c("log","H"),levels=c("data","log","H")),each=nrow(results.df)/2))
colnames(results.df)[5] = "score_type"
#--------------------------------------------------------------------------------------------
#### Artificially put the observations into the same dataframe in order to plot them using facets
temp = subset(results.df,Model==1&repl==1)
temp$score = c(observations[1,],observations[2,])
temp$Model = factor(0)
temp$repl = factor(rep(c(1,2),each=41))
temp$score_type = factor("data")
results.df = rbind(results.df,temp)
#--------------------------------------------------------------------------------------------
#### Artificially arrange observations in the dataframe in order to plot SEGMENTS within the facets
results.df$yend = 0
results.df = rbind(results.df, data.frame(score = observations[1,],
                                          time = 1:41,
                                          Model = factor(0),
                                          repl = factor(0),
                                          score_type = factor("data"),
                                          yend = observations[2,]))
#--------------------------------------------------------------------------------------------
case_label <- list(
  'data'=expression(paste("Observations",sep="")),
  'log'=expression(paste("log-score / Time",sep="")),
  'H'=expression(paste("H-score",sep=""))
)
case_labeller <- function(variable,value){
  return(case_label[value])
}
#------------------------------------------------------------------------------
labels.df = data.frame(x = rep(43.5,6), y = c(14.1,13.6,13.1,-3.55e-3,-4.1e-3,-4.65e-3),
                       text = rep(c("Model 1","Model 2","Model 3"),2),
                       type = factor(rep(c("Model 1","Model 2","Model 3"),2)),
                       score_type = factor(rep(c("log","H"),each=3)))
hlines.df = data.frame(yintercept = 12.8,
                       score_type = factor("log",levels = c("data","log","H")))
colors = wes_palette("Darjeeling")[c(1,3,5)]
axis_titlesize = 22
axis_ticktextsize = 15
ggplot() +
  xlab("Time (number of observations)") +
  ylab("") +
  xlim(0,43.5) +
  scale_fill_manual(values = colors) +
  facet_grid(score_type~., scales="free", labeller = case_labeller) +
  # geom_line(data = subset(results.df,score_type=="data"), aes(time,score,group=repl), color = wes_palette("Darjeeling")[1]) +
  geom_segment(data = subset(results.df,score_type=="data"&repl==0),aes(x = time, y = score, xend = time, yend = yend),linetype="dotted",colour="black") +
  geom_point(data = subset(results.df,score_type=="data"), aes(time,score), color = wes_palette("Darjeeling2")[2]) +
  geom_line(data = subset(results.df,score_type=="H"),aes(time, score, color = Model, group=interaction(Model,repl)),alpha=0.7) +
  geom_line(data = subset(results.df,score_type=="log"),aes(time, score/time, color = Model, group=interaction(Model,repl)),alpha=0.7) +
  geom_hline(data = hlines.df, aes(yintercept = yintercept), alpha=0) +
  stat_summary(data = subset(results.df,score_type=="H"),aes(time, score, color = Model,group=interaction(Model),shape = Model),geom="point", fun.y=mean,size=2.5) +
  stat_summary(data = subset(results.df,score_type=="H"),aes(time, score, color = Model,group=interaction(Model)),geom="line", fun.y=mean, size = 1) +
  stat_summary(data = subset(results.df,score_type=="log"),aes(time, score/time, color = Model,group=interaction(Model),shape = Model),geom="point", fun.y=mean,size=2.5) +
  stat_summary(data = subset(results.df,score_type=="log"),aes(time, score/time, color = Model,group=interaction(Model)),geom="line", fun.y=mean, size = 1) +
  scale_color_manual(values = rep(colors,2)) +
  geom_label(data = labels.df, aes(x,y,label = text,color=rep(colors[3:1],2),fontface="bold"),size=5) +
  theme(axis.text.x = element_text(size = axis_ticktextsize),
        axis.text.y = element_text(size = axis_ticktextsize),
        axis.title.x = element_text(size = axis_titlesize, margin=margin(20,0,0,0)),
        axis.title.y = element_text(size = axis_titlesize, angle = 90, margin = margin(0,20,0,0)),
        strip.text.y = element_text(size = axis_titlesize, colour = "black"),
        strip.background = element_rect(fill="gray88"),
        panel.background = element_rect(fill="gray95",linetype = "solid", colour="white"),
        legend.position = "none")
# ggsave("example_kangaroos_15_by_9.png",dpi = 300,width = 15,height = 9)
# ggsave("example_kangaroos_15_by_9.pdf",dpi = 300,width = 15,height = 9)
