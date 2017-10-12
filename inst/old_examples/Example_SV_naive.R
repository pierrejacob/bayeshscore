##################################################################################################
# Example 2: Levy-driven Stochastic Volatility models (Chopin, Jacob, Papaspiliopoulos, 2013)
# illustration of DDE vs filtering means
##################################################################################################
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(wesanderson)
set.seed(7)
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$Nx = 2^7
algorithmic_parameters$Nx_max = 2^10
algorithmic_parameters$nmoves = 1
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$use_dde = TRUE
algorithmic_parameters$dde_options = list(Ny = 10^4,
                                          sigma2_order0 = 0.001,
                                          sigma2_order1 = 0.002,
                                          sigma2_order2 = 0.01,
                                          nb_steps = Inf)
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Case 1: simulated data from model 1 (single factor Levy-driven)
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# simulate observations
nobservations = 1
timesteps = 1:nobservations
theta = c(0, 0, 0.5, 0.0625, 0.01)
observations = simulateData(get_model_SVLevy_singlefactor(timesteps),theta,nobservations)$Y
#--------------------------------------------------------------------------------------------
# define models
model = function(i){
  if(i==1){return(get_model_SVLevy_singlefactor(timesteps))}
  if(i==2){return(get_model_SVLevy_multifactor_noleverage(timesteps))}
  if(i==3){return(get_model_SVLevy_multifactor_withleverage(timesteps))}
}
#--------------------------------------------------------------------------------------------
repl = 100 #number of replications
# registerDoParallel(cores=5) #number of workers in parallel
models_to_run = c(1)
#--------------------------------------------------------------------------------------------
### Compute logevidence and hscore for each model
results_all = data.frame()
post_all = list()
#--------------------------------------------------------------------------------------------
# # Monitor progress in parallel via log file
# logfilename = "results.log"
# writeLines(c(""), logfilename)
# sink(logfilename, append = TRUE)
# #--------------------------------------------------------------------------------------------
# for (m in models_to_run){
#   cat("Model ",toString(m)," started at:", toString(Sys.time()))
#   results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
#     sink(logfilename, append = TRUE) # Monitor progress in parallel via log file
#     hscore(observations, model(m), algorithmic_parameters)
#   }
#   post = data.frame()
#   for (r in 1:repl){
#     results_all = rbind(results_all,data.frame(logevidence = results[[r]]$logevidence,
#                                                hscore = results[[r]]$hscore,
#                                                hscoreDDE = results[[r]]$hscoreDDE,
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
results = list()
for (m in models_to_run){
  cat("Model ",toString(m)," started at:", toString(Sys.time()))
  for (i in 1:repl){
    print(paste("Replication",toString(i)))
    results[[i]] = hscore(observations, model(m), algorithmic_parameters)
  }
  post = data.frame()
  for (r in 1:repl){
    results_all = rbind(results_all,data.frame(logevidence = results[[r]]$logevidence,
                                               hscore = results[[r]]$hscore,
                                               hscoreDDE = results[[r]]$hscoreDDE,
                                               time = 1:nobservations,
                                               model = factor(m),
                                               repl = r))
    post = rbind(post,data.frame(t(results[[r]]$thetas),
                                 W = results[[r]]$normw,
                                 repl = r))
  }
  post_all[[m]] = post
}
# #--------------------------------------------------------------------------------------------
# # close log file
# sink()
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Generate plots
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
# Check the log-evidence
ggplot(results_all) +
  geom_boxplot(aes(1,logevidence))
# ggsave("example_6_SV_logevidence.png",dpi = 300,width = 10,height = 5)
#--------------------------------------------------------------------------------------------
# Check the hyvarinen score
htype = factor(c("Filtering","DDE"),levels = c("Filtering","DDE"))
hscore.df = data.frame(type = rep(htype,each=repl))

# Define y-axis transformation (signed square root) for readability
signed_sqrt <- scales::trans_new("signed_sqrt", transform = function(x) sign(x)*sqrt(abs(x)),
                                 inverse = function(x) sign(x)*(abs(x)^2))
colors = wes_palette("FantasticFox")[c(5,3)]
#-------------------------
labels.df = data.frame(x = htype, y = rep(-5*10^5,2),
                       text = c("Filtering","DDE"),
                       type = htype)
# Hyvarinen factor
# seed = sample(1:1000,1)
# print(seed)
set.seed(265)
# set.seed(seed)
hscore.df$hscore = c(sample(results_all$hscore),sample(results_all$hscoreDDE))
ggplot() +
  geom_jitter(data=hscore.df,aes(type, hscore,group=type,color=type,shape=type), alpha= 0.8,size = 2,width = 1, height = 0) +
  ylab("Prequenial Hyvärinen score") +
  xlab("Estimation method") +
  geom_label(data = labels.df, aes(x,y,label = text,color=type), color = colors, fontface = "bold") +
  theme(legend.position="none") +
  # scale_y_continuous(trans = signed_log) +
  scale_y_continuous(trans = signed_sqrt) +
  theme(legend.position="none") +
  scale_color_manual(values = colors) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave(paste("example_SV_hscore_naive_",toString(seed),".png",sep=""),dpi = 300,width = 10,height = 5)


