##################################################################################################
# Example: Levy-driven Stochastic Volatility models (Chopin, Jacob, Papaspiliopoulos, 2013)
# (see also Barndorff-Nielsen and Shephard, 2001)
##################################################################################################
# rm(list = ls())
library(bayeshscore)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(wesanderson)
set.seed(19)
#--------------------------------------------------------------------------------------------
# Monitor progress in parallel via log file
# setwd("C:/Users/shao/Desktop/HyvarinenSSM")
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
algorithmic_parameters$nmoves = 10
algorithmic_parameters$verbose = TRUE
# use direct density-estimation
algorithmic_parameters$use_dde = TRUE
algorithmic_parameters$dde_options = list(Ny = 2^10, nb_steps = Inf,
                                          sigma2_order0 = 0.1, sigma2_order1 = 0.1, sigma2_order2 = 0.1)
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
model = function(i){
  if(i==1){return(get_model_SVLevy_singlefactor(timesteps))}
  if(i==2){return(get_model_SVLevy_multifactor_noleverage(timesteps))}
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
observations = readRDS("../Harvard/MyHarvard/_Research/Model Selection for SSM (Pierre and Jie)/Simulation/Example 6 - Stochastic Volatility (Chopin et al., 2013)/T=1000_1/observations_1000.rds")

# Plot observations
ggplot() +
  geom_line(aes(1:ncol(observations),observations[1,]), color = "forestgreen") +
  ylab("Observations") +
  scale_x_continuous(breaks=seq(0,1000,200)) +
  xlab("Time")
# ggsave("example_SV_observations.png",width = 10, height = 5,dpi = 300)

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
  # If error due to lack of RAM, re-run until completion (hacky but the while loop only kicks in when R runs out of RAM)
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
my_theme = theme(axis.text.x = element_text(size = axis_ticktextsize),
                 axis.text.y = element_text(size = axis_ticktextsize),
                 axis.title.x = element_text(size = axis_titlesize, margin=margin(20,0,0,0)),
                 axis.title.y = element_text(size = axis_titlesize, angle = 90, margin = margin(0,20,0,0)),
                 strip.text.y = element_text(size = axis_titlesize, colour = "black"),
                 strip.text.x = element_text(size = axis_titlesize, colour = "black",margin = margin(0.3,0,0.3,0, "cm")),
                 strip.background = element_rect(fill="gray88"),
                 panel.background = element_rect(fill="gray95",linetype = "solid", colour="white"),
                 legend.position = "none")
# Check posterior for each model
post_plot_all = list()
axis_titlesize = 24
axis_ticktextsize = 15
xlabels = c(expression(mu), expression(beta), expression(xi), expression(omega^2),
            expression(lambda[1]), expression(lambda[2]), "w", expression(rho[1]), expression(rho[2]))
colors = c(wes_palette("Darjeeling2")[c(2,3)],wes_palette("Darjeeling")[2])
# Generate plots
for (m in models_to_run){
  dimtheta = model(m)$dimtheta
  post = data.frame(thetas = unlist(post_all[[m]][,1:dimtheta]),
                    repl = rep(post_all[[m]]$repl,dimtheta),
                    w = rep(post_all[[m]]$W,dimtheta),
                    variable = rep(factor(1:dimtheta),each=length(post_all[[m]]$repl)),
                    type = factor("Posterior"))
  if (m==1){levels(post$variable) = xlabels[1:5]; nbcol = 5}
  if (m==2){levels(post$variable) = xlabels[1:7]; nbcol = 7}

  local({m = m;
  post_plot_all[[m]] <<-ggplot(post) +
    geom_density(aes(thetas, weight = w, group = interaction(variable,repl)),adjust = 2,col=colors[m],size=1,alpha=0.3) +
    facet_wrap( ~ variable, ncol=nbcol, scales="free",labeller=label_parsed) + xlab("") + ylab("") +
    my_theme})
}
post_plot_all[[1]]
post_plot_all[[2]]
# ggsave("example_SV_posterior_model_1_24_by_6.png",plot = post_plot_all[[1]],dpi = 300,width = 24,height = 6)
# ggsave("example_SV_posterior_model_2_24_by_6.png",plot = post_plot_all[[2]],dpi = 300,width = 24,height = 6)



#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------

a = 0.95
# Check the log-evidence
ggplot(results_all,aes(time,  -logevidence, color = model)) +
  ylab("- log evidence") +
  # guides(shape = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_x_continuous(breaks=seq(0,1000,200)) +
  geom_line(aes(group=interaction(model,repl)),linetype="dashed",alpha=0.8) +
  # stat_summary(aes(group=model,fill=model),geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=a),alpha=0.3) +
  stat_summary(aes(group=model,shape = model),geom="point", fun.y=mean,size=2) +
  stat_summary(aes(group=model),geom="line", fun.y=mean, size = 1)
# ggsave("example_SV_logevidence.png",dpi = 300,width = 10, height = 5)
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Check the h-score DDE
ggplot(subset(results_all),aes(time, hscoreDDE/time, color = model)) +
  ylab("Hyvarinen score / time") +
  # guides(shape = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_x_continuous(breaks=seq(0,1000,200)) +
  geom_line(aes(group=interaction(model,repl)),linetype="dashed",alpha=0.5) +
  # stat_summary(aes(group=model,fill=model),geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=a),alpha=0.3) +
  stat_summary(aes(group=model,shape = model),geom="point", fun.y=mean,size=2) +
  stat_summary(aes(group=model),geom="line", fun.y=mean, size = 1)
# ggsave("example_SV_hscore_rescaled.png",width = 10, height = 5,dpi = 300)
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
colors = wes_palette("Darjeeling2")[c(2,3)]
labels.df = data.frame(x = rep(1075,2), y = c(-2900,-2680),
                       text = c("Model 1","Model 2"),
                       type = factor(c("Model 1","Model 2")))
x_resolution = (1:nobservations)[(1:nobservations)%%10==0]
# x_resolution = (1:nobservations)
# Check the h-score DDE
ggplot(subset(results_all,time%in%x_resolution),aes(time, hscoreDDE, color = model)) +
  ylab("Prequential Hyvärinen score") +
  xlab("Time") +
  xlim(0,1100) +
  scale_x_continuous(breaks=seq(0,1000,200)) +
  geom_label(data = labels.df, aes(x,y,label = text,color=type), color = colors, fontface = "bold") +
  theme(legend.position="none") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  geom_line(aes(group=interaction(model,repl)),linetype="solid",alpha=0.5) +
  # stat_summary(aes(group=model,fill=model),geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=a),alpha=0.3) +
  stat_summary(aes(group=model,shape = model),geom="point", fun.y=mean,size=2) +
  stat_summary(aes(group=model),geom="line", fun.y=mean, size = 1)
# ggsave("example_SV_hscore_every10steps.png",width = 10, height = 5,dpi = 300)
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
x_resolution = (1:nobservations)[(1:nobservations)%%10==0]
# x_resolution = (1:nobservations)
centered = results_all
centered$hscoreDDE = centered$hscoreDDE - rep(rowMeans(sapply(1:repl,function(i)subset(centered,model==1&repl==i)$hscoreDDE)),repl*length(unique(centered$model)))
labels.df = data.frame(x = rep(1075,2), y = c(0,mean(subset(centered,time==1000 & model==2)$hscoreDDE)),
                       text = c("Model 1","Model 2"),
                       type = factor(c("Model 1","Model 2")))
# Check the h-score DDE
ggplot(subset(centered,time%in%x_resolution),aes(time, hscoreDDE, color = model)) +
  ylab("Prequential Hyvärinen score (centered)") +
  # guides(shape = guide_legend(override.aes = list(size=2))) +
  xlim(0,1100) +
  xlab("Time") +
  geom_label(data = labels.df, aes(x,y,label = text,color=type), color = colors, fontface = "bold") +
  theme(legend.position="none") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  geom_line(aes(group=interaction(model,repl)),linetype="solid",alpha=0.5) +
  # stat_summary(aes(group=model,fill=model),geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=a),alpha=0.3) +
  stat_summary(aes(group=model,shape = model),geom="point", fun.y=mean,size=2) +
  stat_summary(aes(group=model),geom="line", fun.y=mean, size = 1) +
  scale_x_continuous(breaks=seq(0,1000,200))
# ggsave("example_SV_hscore_shift_every10steps.png",width = 10, height = 5,dpi = 300)
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#### Artificially put all the things to be plotted into a single dataframe in order to use facets
BHfactors.df = data.frame()
for (r in 1:repl){
  BHfactors.df = rbind(BHfactors.df, data.frame(value = subset(centered,model==2&repl==r)$hscoreDDE-subset(centered,model==1&repl==r)$hscoreDDE,
                                                time = 1:nobservations,
                                                repl = r,
                                                type = factor("H",levels = c("data","log","H"))))
  BHfactors.df = rbind(BHfactors.df, data.frame(value = -subset(centered,model==2&repl==r)$logevidence+subset(centered,model==1&repl==r)$logevidence,
                                                time = 1:nobservations,
                                                repl = r,
                                                type = factor("log",levels = c("data","log","H"))))
}
BHfactors.df = rbind(BHfactors.df, data.frame(value = c(observations),
                                              time = 1:nobservations,
                                              repl = 0,
                                              type = factor("data",levels = c("data","log","H"))))
hlines.df = data.frame(yintercept = rep(0,2),
                       type = factor(c("log","H"),levels = c("data","log","H")))
#--------------------------------------------------------------------------------------------
case_label <- list(
  'data'=expression(paste("Observations",sep="")),
  'log'=expression(paste("log-BF 1 vs. 2",sep="")),
  'H'=expression(paste("HF 1 vs. 2",sep=""))
)
case_labeller <- function(variable,value){
  return(case_label[value])
}
colors = c("forestgreen","tomato",wes_palette("Darjeeling2")[2])
mean_size = 1.5
axis_titlesize = 18
axis_ticktextsize = 15
#--------------------------------------------------------------------------------------------
ggplot() +
  xlab("Time (number of observations)") +
  ylab("") +
  # scale_fill_manual(values = colors) +
  facet_grid(type~., scales="free", labeller = case_labeller) +
  # geom_line(data = subset(results.df,score_type=="data"), aes(time,score,group=repl), color = wes_palette("Darjeeling")[1]) +
  geom_line(data = subset(BHfactors.df,type=="data"),aes(time, value),alpha=1,color = colors[1]) +
  geom_hline(data = hlines.df, aes(yintercept = yintercept), linetype ="dashed") +
  geom_line(data = subset(BHfactors.df,type=="H"),aes(time, value, group=repl),alpha=0.7,color = colors[3]) +
  geom_line(data = subset(BHfactors.df,type=="log"),aes(time, value, group=repl),alpha=0.7,color = colors[2]) +
  stat_summary(data = subset(BHfactors.df,type=="H"),aes(time, value) ,geom="line", fun.y=mean, size = mean_size,color = colors[3]) +
  stat_summary(data = subset(BHfactors.df,type=="log"),aes(time, value) ,geom="line", fun.y=mean, size = mean_size,color = colors[2]) +
  my_theme
# ggsave("example_SV_12_by_9.png",dpi = 300,width = 12,height = 9)
# ggsave("example_SV_15_by_9.png",dpi = 300,width = 15,height = 9)



# #### Plot prior (WARNING: these are specific to a particular choice of prior, xmin and xmax were set manually
# #### for better readability of the plots)
# mu0mu = 0
# sigma02mu = 10
# mu0beta = 0
# sigma02beta = 10
# r0xi = 1/5
# r0w2 = 1/5
# r0lambda = 1
# ### Model 1
# M = 5000
# axis_titlesize = 24
# axis_ticktextsize = 15
# my_theme = theme(axis.text.x = element_text(size = axis_ticktextsize),
#                  axis.text.y = element_text(size = axis_ticktextsize),
#                  axis.title.x = element_text(size = axis_titlesize, margin=margin(20,0,0,0)),
#                  axis.title.y = element_text(size = axis_titlesize, angle = 90, margin = margin(0,20,0,0)),
#                  strip.text.y = element_text(size = axis_titlesize, colour = "black"),
#                  strip.text.x = element_text(size = axis_titlesize, colour = "black",margin = margin(0.3,0,0.3,0, "cm")),
#                  strip.background = element_rect(fill="gray88"),
#                  panel.background = element_rect(fill="gray95",linetype = "solid", colour="white"),
#                  legend.position = "none")
# xmin = -0.15; xmax = 0.17
# prior1 = data.frame(thetas = seq(xmin,xmax,length.out = M),
#                     dens = dnorm(seq(xmin,xmax,length.out = M),mu0mu, sqrt(sigma02mu)),
#                     variable = factor(1),
#                     type = factor("Prior"))
# xmin = -0.5; xmax = 0.45
# prior1 = rbind(prior1, data.frame(thetas = seq(xmin,xmax,0.0005),
#                                   dens = dnorm(seq(xmin,xmax,0.0005),mu0beta, sqrt(sigma02beta)),
#                                   variable = factor(2),
#                                   type = factor("Prior")))
# xmin = 0; xmax = 2.8
# prior1 = rbind(prior1, data.frame(thetas = seq(xmin,xmax,length.out = M),
#                                   dens = dexp(seq(xmin,xmax,length.out = M),r0xi),
#                                   variable = factor(3),
#                                   type = factor("Prior")))
# xmin = 0; xmax = 6
# prior1 = rbind(prior1, data.frame(thetas = seq(xmin,xmax,length.out = M),
#                                   dens = dexp(seq(xmin,xmax,length.out = M),r0w2),
#                                   variable = factor(4),
#                                   type = factor("Prior")))
# xmin = 0; xmax = 0.04
# prior1 = rbind(prior1, data.frame(thetas = seq(xmin,xmax,length.out = M),
#                                   dens = dexp(seq(xmin,xmax,length.out = M),r0lambda),
#                                   variable = factor(5),
#                                   type = factor("Prior")))
# levels(prior1$variable) = xlabels[1:5]
# ylim.df = data.frame(y = c(0.12,0.13,0.1,0.15,0,1,0,1,0,2),
#                      variable = factor(rep(1:5,each=2)))
# levels(ylim.df$variable) = xlabels[1:5]
# ggplot(prior1) +
#   geom_line(aes(thetas, dens), col=colors[1],size=1,alpha=0.8) +
#   facet_wrap( ~ variable, ncol=5, scales="free",labeller=label_parsed) + xlab("") + ylab("") +
#   geom_hline(data = ylim.df, aes(yintercept = y), alpha = 0) +
#   my_theme
# # ggsave("example_SV_prior_model_1_24_by_6.png",dpi = 300,width = 24,height = 6)
#
#
# r0lambda1 = 1
# r0lambda2_1 = 1/2
# alpha0w = 1
# beta0w = 1
# ### Model 2
# M = 5000
# xmin = -0.2; xmax = 0.2
# prior2 = data.frame(thetas = seq(xmin,xmax,length.out = M),
#                     dens = dnorm(seq(xmin,xmax,length.out = M),mu0mu, sqrt(sigma02mu)),
#                     variable = factor(1),
#                     type = factor("Prior"))
# xmin = -0.5; xmax = 0.5
# prior2 = rbind(prior2, data.frame(thetas = seq(xmin,xmax,0.0005),
#                                   dens = dnorm(seq(xmin,xmax,0.0005),mu0beta, sqrt(sigma02beta)),
#                                   variable = factor(2),
#                                   type = factor("Prior")))
# xmin = 0.2; xmax = 0.8
# prior2 = rbind(prior2, data.frame(thetas = seq(xmin,xmax,length.out = M),
#                                   dens = dexp(seq(xmin,xmax,length.out = M),r0xi),
#                                   variable = factor(3),
#                                   type = factor("Prior")))
# xmin = 0; xmax = 0.8
# prior2 = rbind(prior2, data.frame(thetas = seq(xmin,xmax,length.out = M),
#                                   dens = dexp(seq(xmin,xmax,length.out = M),r0w2),
#                                   variable = factor(4),
#                                   type = factor("Prior")))
# xmin = 0; xmax = 0.05
# prior2 = rbind(prior2, data.frame(thetas = seq(xmin,xmax,length.out = M),
#                                   dens = dexp(seq(xmin,xmax,length.out = M),r0lambda1),
#                                   variable = factor(5),
#                                   type = factor("Prior")))
# xmin = 0; xmax = 14
# draws = model(2)$rprior(10^6)[6,]
# # ggplot() + geom_histogram(aes(x = draws, y = ..density..))
# kde = density(draws,bw = 1.5)
# prior2 = rbind(prior2, data.frame(thetas = kde$x[kde$x>xmin & kde$x<xmax],
#                                   dens = kde$y[kde$x>xmin & kde$x<xmax],
#                                   variable = factor(6),
#                                   type = factor("Prior")))
# xmin = 0.35; xmax = 1
# prior2 = rbind(prior2, data.frame(thetas = seq(xmin,xmax,length.out = M),
#                                   dens = dbeta(seq(xmin,xmax,length.out = M), alpha0w, beta0w),
#                                   variable = factor(7),
#                                   type = factor("Prior")))
# levels(prior2$variable) = xlabels[1:7]
# ylim.df = data.frame(y = c(0.12,0.13,0.1,0.15,0,1,0,1,0,2,0,0.2,0,1.5),
#                      variable = factor(rep(1:7,each=2)))
# levels(ylim.df$variable) = xlabels[1:7]
# ggplot(prior2) +
#   geom_line(aes(thetas, dens), col=colors[2],size=1,alpha=0.8) +
#   facet_wrap( ~ variable, ncol=7, scales="free",labeller=label_parsed) + xlab("") + ylab("") +
#   geom_hline(data = ylim.df, aes(yintercept = y), alpha = 0) +
#   my_theme
# # ggsave("example_SV_prior_model_2_24_by_6.png",dpi = 300,width = 24,height = 6)
