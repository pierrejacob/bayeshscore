##################################################################################################
# Example - AR(1) vs MA(1)
##################################################################################################
library(bayeshscore)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(wesanderson)
set.seed(19)

##########################################################################
# WARNING: Models 1 and 2 here correspond to Models 1 and 3 in the paper.
##########################################################################

# Define model
nu0 = 1
sigma02 = 1
nb_models = 2
model = function(i){
  if (i==1){return(get_model_AR1(nu0, sigma02))}
  if (i==2){return(get_model_ARMA(0,1,nu0, sigma02))}
}
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$Nx = 2^7
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$nmoves = 1
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------
# Monitor progress in parallel via log file
# setwd("C:/Users/shao/Desktop/bayeshscore")
logfilename = "results.log"
writeLines(c(""), logfilename)
sink(logfilename, append = TRUE)
#--------------------------------------------------------------------------------------------
nobservations = 1000
##################################################################################################
# Case 1: true model = AR(1)
##################################################################################################
true_model = 1
true_theta = c(0.6,0.8^2)
observations1 = simulateData(model(true_model),true_theta,nobservations)
#--------------------------------------------------------------------------------------------
results_all1 = data.frame()
### Compute logevidence and hscore for each model
for (m in 1:nb_models){
  results = foreach(i=1:repl,.packages=c('bayeshscore'),.verbose = TRUE) %dorng% {
    sink(logfilename, append = TRUE) # Monitor progress in parallel via log file
    hscore(observations1, model(m), algorithmic_parameters)
  }
  for (r in 1:repl){
    results_all1 = rbind(results_all1,data.frame(logevidence = results[[r]]$logevidence,
                                                 hscore = results[[r]]$hscore,
                                                 time = 1:nobservations,
                                                 model = factor(m),
                                                 repl = r))
  }
}
##################################################################################################
# Case 2: true model = MA(1)
##################################################################################################
true_model = 2
true_theta = c(0.5,1)
observations2 = simulateData(model(true_model),true_theta,nobservations)$Y
#--------------------------------------------------------------------------------------------
results_all2 = data.frame()
### Compute logevidence and hscore for each model
for (m in 1:nb_models){
  results = foreach(i=1:repl,.packages=c('bayeshscore'),.verbose = TRUE) %dorng% {
    sink(logfilename, append = TRUE) # Monitor progress in parallel via log file
    hscore(observations2, model(m), algorithmic_parameters)
  }
  for (r in 1:repl){
    results_all2 = rbind(results_all2,data.frame(logevidence = results[[r]]$logevidence,
                                                 hscore = results[[r]]$hscore,
                                                 time = 1:nobservations,
                                                 model = factor(m),
                                                 repl = r))
  }
}
#--------------------------------------------------------------------------------------------
# close log file
sink()
##################################################################################################
# Compare model 1 and model 3 in each case
##################################################################################################
results_all = list(results_all1,results_all2)
logbayesfactors = data.frame()
h_factors = data.frame()
BF_plots = list()
HF_plots = list()
for (r in 1:repl) {
  for (i in 1:nb_models) {
    results = results_all[[i]]
    logbayes_factor = subset(results,model==1&repl==r)$logevidence - subset(results,model==2&repl==r)$logevidence
    logbayesfactors = rbind(logbayesfactors,data.frame(case = factor(i),
                                                       time = 1:nobservations,
                                                       repl = r,
                                                       logbayesfactor = logbayes_factor,
                                                       type = factor(paste("Case ",toString(i)))))
    h_factor = subset(results,model==2&repl==r)$hscore - subset(results,model==1&repl==r)$hscore
    h_factors = rbind(h_factors,data.frame(case = factor(i),
                                           time = 1:nobservations,
                                           repl = r,
                                           hfactor = h_factor,
                                           type = factor(paste("Case ",toString(i)))))
  }
}
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Generate plots for paper
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Bayes factor
ggplot(logbayesfactors) +
  geom_line(aes(time, logbayesfactor, color = case, group = interaction(case,repl))) +
  geom_hline(yintercept = 0,linetype="dotted",size=1) +
  ylab("log Bayes factor  [1 vs 3]") + facet_grid(. ~ type) + xlab("Number of observations") +
  # guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave("example_AR1_MA1_log_BF_1_vs_2.png",dpi = 300,width = 10,height = 5)

case_label <- list(
  'Case 1'=expression(paste("Case 1: ",M[1]," is well-specified",sep="")),
  'Case 2'=expression(paste("Case 2: ",M[3]," is well-specified",sep=""))
)
case_labeller <- function(variable,value){
  return(case_label[value])
}
# Hyvarinen factor
ggplot(h_factors) +
  geom_line(aes(time, hfactor, color = case, group = interaction(case,repl))) +
  geom_hline(yintercept = 0,linetype="dotted",size=1) +
  xlab("Number of observations") +
  ylab("H-factor  [1 vs 3]") +
  facet_grid(. ~ type, labeller = case_labeller) +
  # guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave("example_AR1_MA1_Hyvarinen_factor_1_vs_2.png",dpi = 300,width = 10,height = 5)

# Hyvarinen factor
labels10by5.df = data.frame(x = rep(250,2), y = c(-10,22),
                       text = c("Case 2","Case 1"),
                       type = factor(c("Case 2","Case 1")))
colors = c(wes_palette("Darjeeling")[c(4,2)])
ggplot() +
  geom_label(data = labels10by5.df, aes(x,y,label = text,color=type), color = colors, fontface = "bold") +
  scale_color_manual(values = colors[2:1]) +
  geom_line(data = h_factors, aes(time, hfactor, color = case, group = interaction(case,repl)),alpha=0.6) +
  geom_hline(yintercept = 0,linetype=2) +
  xlab("Number of observations") +
  ylab("H-factor  [1 vs. 3]") +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave("example_AR1_MA1_Hyvarinen_factor_1_vs_2_10_by_5.png",dpi = 300,width = 10,height = 5)

# Hyvarinen factor
labels5by5.df = data.frame(x = rep(187,2), y = c(-15,22),
                       text = c("Case 2","Case 1"),
                       type = factor(c("Case 2","Case 1")))
ggplot() +
  geom_label(data = labels5by5.df, aes(x,y,label = text,color=type), color = colors, fontface = "bold") +
  scale_color_manual(values = colors[2:1]) +
  geom_line(data = h_factors, aes(time, hfactor, color = case, group = interaction(case,repl)),alpha=0.6) +
  geom_hline(yintercept = 0,linetype=2) +
  xlab("Number of observations") +
  ylab("H-factor  [1 vs. 3]") +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(axis.text.x = element_text(size = axis_ticktextsize),
        axis.text.y = element_text(size = axis_ticktextsize),
        axis.title.x = element_text(size = axis_titlesize, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = axis_titlesize, angle = 90, margin = margin(0,10,0,0)),
        strip.text.x = element_text(size = axis_titlesize, colour = "black"),
        strip.background = element_rect(fill="gray88"),
        panel.background = element_rect(fill="gray95",linetype = "solid", colour="white"),
        legend.position = "none")
# ggsave("example_AR1_MA1_Hyvarinen_factor_1_vs_2_5_by_5.png",dpi = 300,width = 5,height = 5)



# log Bayes factor
labels10by5.df = data.frame(x = rep(250,2), y = c(-8,12),
                            text = c("Case 2","Case 1"),
                            type = factor(c("Case 2","Case 1")))
colors = c(wes_palette("Darjeeling")[c(4,2)])
ggplot() +
  geom_label(data = labels10by5.df, aes(x,y,label = text,color=type), color = colors, fontface = "bold") +
  scale_color_manual(values = colors[2:1]) +
  geom_line(data = logbayesfactors, aes(time, logbayesfactor, color = case, group = interaction(case,repl)),alpha=0.6) +
  geom_hline(yintercept = 0,linetype=2) +
  xlab("Number of observations") +
  ylab("log-Bayes factor  [1 vs. 3]") +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  ylim(c(min(h_factors$hfactor),max(h_factors$hfactor))) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave("example_AR1_MA1_logBF_1_vs_2_10_by_5.png",dpi = 300,width = 10,height = 5)


# log Bayes factor
labels5by5.df = data.frame(x = rep(187,2), y = c(-10,13),
                           text = c("Case 2","Case 1"),
                           type = factor(c("Case 2","Case 1")))
ggplot() +
  geom_label(data = labels5by5.df, aes(x,y,label = text,color=type), color = colors, fontface = "bold") +
  scale_color_manual(values = colors[2:1]) +
  geom_line(data = logbayesfactors, aes(time, logbayesfactor, color = case, group = interaction(case,repl)),alpha=0.6) +
  geom_hline(yintercept = 0,linetype=2) +
  xlab("Number of observations") +
  ylab("log-Bayes factor  [1 vs. 3]") +
  ylim(c(min(h_factors$hfactor),max(h_factors$hfactor))) +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(axis.text.x = element_text(size = axis_ticktextsize),
        axis.text.y = element_text(size = axis_ticktextsize),
        axis.title.x = element_text(size = axis_titlesize, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = axis_titlesize, angle = 90, margin = margin(0,10,0,0)),
        strip.text.x = element_text(size = axis_titlesize, colour = "black"),
        strip.background = element_rect(fill="gray88"),
        panel.background = element_rect(fill="gray95",linetype = "solid", colour="white"),
        legend.position = "none")
# ggsave("example_AR1_MA1_logBF_1_vs_2_5_by_5.png",dpi = 300,width = 5,height = 5)



#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
crit.df = data.frame()
crit.df = rbind(crit.df,data.frame(logbayesfactors[-4], value = logbayesfactors$logbayesfactor, crit = factor("LBF")))
crit.df = rbind(crit.df,data.frame(h_factors[-4], value = h_factors$hfactor, crit = factor("HF")))
relabelling = function(x){
  if(x==levels(crit.df$type[1])[1]){return (factor(1))}
  if(x==levels(crit.df$type[1])[2]){return (factor(2))}
}
crit.df$type = sapply(crit.df$type, function(x)relabelling(x))
#--------------------------------------------------------------------------------------------
case_label <- list(
  '1'=expression(paste("Case 1: non-nested, ",M[1]," well-specified, ", M[3]," misspecified",sep="")),
  '2'=expression(paste("Case 2: non-nested, ",M[1]," misspecified, ", M[3]," well-specified",sep=""))
)
case_labeller <- function(variable,value){
  return(case_label[value])
}
labels.df = data.frame(x = c(rep(900,2),rep(900,2)),
                       y = c(125,32,-40,-20),
                       text = rep(c("HF 1 vs. 3","log-BF 1 vs. 3"),4),
                       crit = factor(rep(c("HF","LBF"),4)),
                       type = factor(rep(1:2,each=2)))
# hline.df = data.frame(y = c(-5,10), type = factor(3))
axis_titlesize = 18
axis_ticktextsize = 15
ggplot() +
  geom_label(data = labels.df, aes(x,y,label = text,color=crit), fontface = "bold") +
  scale_color_manual(values = colors[2:1]) +
  geom_line(data = crit.df, aes(time, value, color = crit, group = interaction(crit,repl)),alpha=0.8) +
  geom_hline(yintercept = 0,linetype=2) +
  # geom_hline(data = hline.df, aes(yintercept = y),linetype=2, alpha = 0) +
  xlab("Number of observations") +
  ylab("") +
  # ylim(c(min(h_factors$hfactor),30+max(h_factors$hfactor))) +
  facet_wrap( ~ type, ncol=2, scales="free", labeller = case_labeller) +
  theme(axis.text.x = element_text(size = axis_ticktextsize),
        axis.text.y = element_text(size = axis_ticktextsize),
        axis.title.x = element_text(size = axis_titlesize, margin=margin(20,0,0,0)),
        axis.title.y = element_text(size = axis_titlesize, angle = 90, margin = margin(0,20,0,0)),
        strip.text.x = element_text(size = axis_titlesize, colour = "black"),
        strip.text.y = element_text(size = axis_titlesize, colour = "black"),
        legend.position = "none",
        strip.background = element_rect(fill="gray88"),
        panel.background = element_rect(fill="gray95",linetype = "solid", colour="white"))
# ggsave("example_AR1_MA1_15_by_6.png",dpi = 300,width = 15,height = 6)

