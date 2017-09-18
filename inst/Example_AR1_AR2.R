##################################################################################################
# Example - AR(1) vs AR(2)
##################################################################################################
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(wesanderson)
set.seed(19)

# Define model
nu0 = 1
sigma02 = 1
nb_models = 2
model = function(i){
  if (i==1){return(get_model_AR1(nu0, sigma02))}
  if (i==2){return(get_model_AR2(nu0, sigma02))}
}

# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$verbose = TRUE
# The remaining algorithmic parameters are set to their default values via util_default.R
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------
nobservations = 1000

#############################################################################################
# Case 1: true model = AR(1)
#############################################################################################
true_model = 1
true_theta = c(0.5,1)
observations1 = simulateData(model(true_model),true_theta,nobservations)
# observations in a matrix of dimensions dimy x nobservations
#--------------------------------------------------------------------------------------------
results_all1 = data.frame()
#--------------------------------------------------------------------------------------------
### Compute logevidence and hscore for each model
for (m in 1:nb_models){
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
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
#############################################################################################
# Case 2: true model = AR(2)
#############################################################################################
true_model = 2
true_theta = c(0.25,0.35,1)
observations2 = simulateData(model(true_model),true_theta,nobservations)
# observations in a matrix of dimensions dimy x nobservations
#--------------------------------------------------------------------------------------------
results_all2 = data.frame()
#--------------------------------------------------------------------------------------------
### Compute logevidence and hscore for each model
for (m in 1:nb_models){
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
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
##############################################################################################
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Generate plots for paper
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Compute the Hyvarinen factor
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

# log Bayes factor
ggplot(logbayesfactors) +
  geom_line(aes(time, logbayesfactor, color = case, group = interaction(case,repl))) +
  geom_hline(yintercept = 0,linetype="dotted",size=1) +
  ylab("log Bayes factor  [1 vs 2]") + facet_grid(. ~ type) + xlab("Number of observations") +
  # guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave("example_AR1_AR2_logBF_1_vs_2.png",dpi = 300,width = 10,height = 5)


# Hyvarinen factor
case_label <- list(
  'Case 1'=expression(paste("Case 1: both are well-specified",sep="")),
  'Case 2'=expression(paste("Case 2: ",M[2]," is well-specified",sep=""))
)
case_labeller <- function(variable,value){
  return(case_label[value])
}
labels.df = data.frame(x = rep(375,2), y = c(-58,18),
                       text = c("Case 4","Case 3"),
                       type = factor(c("Case 4","Case 3")))
colors = c(wes_palette("Darjeeling")[1],"mediumblue")
ggplot() +
  geom_label(data = labels.df, aes(x,y,label = text,color=type), color = colors, fontface = "bold") +
  theme(legend.position="none") +
  scale_color_manual(values = colors[2:1]) +
  geom_line(data = h_factors, aes(time, hfactor, color = case, group = interaction(case,repl)),alpha=0.6) +
  geom_hline(yintercept = 0,linetype=2) +
  xlab("Number of observations") +
  ylab("Hyvärinen factor  [1 vs 2]") +
  ylim(c(min(h_factors$hfactor),20+max(h_factors$hfactor))) +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave("example_AR1_AR2_Hyvarinen_factor_1_vs_2_10_by_5.png",dpi = 300,width = 10,height = 5)
# ggsave("example_AR1_AR2_Hyvarinen_factor_1_vs_2_5_by_5.png",dpi = 300,width = 5,height = 5)


# log Bayes factor
case_label <- list(
  'Case 1'=expression(paste("Case 1: both are well-specified",sep="")),
  'Case 2'=expression(paste("Case 2: ",M[2]," is well-specified",sep=""))
)
case_labeller <- function(variable,value){
  return(case_label[value])
}
labels.df = data.frame(x = rep(375,2), y = c(-35,15),
                       text = c("Case 4","Case 3"),
                       type = factor(c("Case 4","Case 3")))
colors = c(wes_palette("Darjeeling")[1],"mediumblue")
ggplot() +
  geom_label(data = labels.df, aes(x,y,label = text,color=type), color = colors, fontface = "bold") +
  theme(legend.position="none") +
  scale_color_manual(values = colors[2:1]) +
  geom_line(data = logbayesfactors, aes(time, logbayesfactor, color = case, group = interaction(case,repl)),alpha=0.6) +
  geom_hline(yintercept = 0,linetype=2) +
  xlab("Number of observations") +
  ylab("log Bayes factor  [1 vs 2]") +
  ylim(c(min(h_factors$hfactor),20+max(h_factors$hfactor))) +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave("example_AR1_AR2_logBF_1_vs_2_10_by_5.png",dpi = 300,width = 10,height = 5)
# ggsave("example_AR1_AR2_logBF_1_vs_2_5_by_5.png",dpi = 300,width = 5,height = 5)
