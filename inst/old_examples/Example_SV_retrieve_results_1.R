
setwd("C:/Users/jianfushao/Documents/Harvard/MyHarvard/_Research/Model Selection for SSM (Pierre and Jie)/Simulation/Example 6 - Stochastic Volatility (Chopin et al., 2013)/T=100")

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
results = readRDS("SV_model_1_2_3.rds")
results_all = results$results_all
post_all = results$post_all
repl = 5


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


# Check the log-evidence
ggplot(results_all,aes(time,-logevidence, color = model)) +
  ylab("- log evidence") +
  # guides(shape = guide_legend(override.aes = list(size=2))) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  geom_line(aes(group=interaction(model,repl)),linetype="dashed",alpha=0.8) +
  # stat_summary(aes(group=model,fill=model),geom="ribbon", fun.data=mean_cl_normal, fun.args=list(conf.int=a),alpha=0.3) +
  stat_summary(aes(group=model,shape = model),geom="point", fun.y=mean,size=2) +
  stat_summary(aes(group=model),geom="line", fun.y=mean, size = 1)
