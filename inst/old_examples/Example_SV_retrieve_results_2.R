

setwd("C:/Users/jianfushao/Documents/Harvard/MyHarvard/_Research/Model Selection for SSM (Pierre and Jie)/Simulation/Example 6 - Stochastic Volatility (Chopin et al., 2013)/T=1000_1")

### Retrieve results from saved runs (SV models)
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(wesanderson)
nobservations = 1000
timesteps = 1:nobservations
# define models
model = function(i){
  if(i==1){return(get_model_SVLevy_singlefactor(timesteps))}
  if(i==2){return(get_model_SVLevy_multifactor_noleverage(timesteps))}
  if(i==3){return(get_model_SVLevy_multifactor_withleverage(timesteps))}
}
models_to_run = c(1,2)
results_all = data.frame()
post_all = vector("list",length(models_to_run))
repl = 5


for (m in models_to_run) {
  for (r in 1:repl) {
    results = readRDS(paste("SV_run_model_",toString(m),"_repl_",toString(r),".rds",sep=""))
    # format results into dataframes
    results_all = rbind(results_all,data.frame(logevidence = results$logevidence,
                                               hscore = results$hscore,
                                               hscoreDDE = results$hscoreDDE,
                                               time = 1:nobservations,
                                               model = factor(m),
                                               repl = r))
    post_all[[m]] = rbind(post_all[[m]],data.frame(t(results$thetas),
                                                   W = results$normw,
                                                   repl = r))

  }
}


#--------------------------------------------------------------------------------------------
# Check posterior for each model
post_plot_all = list()
xlabels = list(xlab(expression(mu)), xlab(expression(beta)), xlab(expression(xi)),
               xlab(expression(omega^2)), xlab(expression(lambda[1])), xlab(expression(lambda[2])),
               xlab("w"), xlab(expression(rho[1])), xlab(expression(rho[2])))
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
ggplot(subset(results_all),aes(time, hscoreDDE/time, color = model)) +
  ylab("Hyvarinen score / time") +
  # guides(shape = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  geom_line(aes(group=interaction(model,repl)),linetype="dashed",alpha=0.5) +
  # stat_summary(aes(group=model,fill=model),geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=a),alpha=0.3) +
  stat_summary(aes(group=model,shape = model),geom="point", fun.y=mean,size=2) +
  stat_summary(aes(group=model),geom="line", fun.y=mean, size = 1)
# ggsave("example_SV_hscore_rescaled.png",width = 10, height = 5,dpi = 300)


colors = wes_palette("Darjeeling2")[c(2,3)]
labels.df = data.frame(x = rep(1075,2), y = c(-2900,-2650),
                       text = c("Model 1","Model 2"),
                       type = factor(c("Model 1","Model 2")))
x_resolution = (1:nobservations)[(1:nobservations)%%10==0]
# x_resolution = (1:nobservations)
# Check the h-score DDE
ggplot(subset(results_all,time%in%x_resolution),aes(time, hscoreDDE, color = model)) +
  ylab("Prequential Hyvärinen score") +
  xlab("Time") +
  xlim(0,1100) +
  geom_label(data = labels.df, aes(x,y,label = text,color=type), color = colors, fontface = "bold") +
  theme(legend.position="none") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  geom_line(aes(group=interaction(model,repl)),linetype="solid",alpha=0.5) +
  # stat_summary(aes(group=model,fill=model),geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=a),alpha=0.3) +
  stat_summary(aes(group=model,shape = model),geom="point", fun.y=mean,size=2) +
  stat_summary(aes(group=model),geom="line", fun.y=mean, size = 1)
# ggsave("example_SV_hscore_every10steps.png",width = 10, height = 5,dpi = 300)


x_resolution = (1:nobservations)[(1:nobservations)%%10==0]
# x_resolution = (1:nobservations)
centered = results_all
centered$hscoreDDE = centered$hscoreDDE - rep(rowMeans(sapply(1:repl,function(i)subset(centered,model==1&repl==i)$hscoreDDE)),repl*length(unique(centered$model)))
labels.df = data.frame(x = rep(1075,2), y = c(0,90),
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
  stat_summary(aes(group=model),geom="line", fun.y=mean, size = 1)
# ggsave("example_SV_hscore_shift_every10steps.png",width = 10, height = 5,dpi = 300)
