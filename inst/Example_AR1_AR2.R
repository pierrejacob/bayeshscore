##################################################################################################
# Example - AR(1) vs AR(2)
##################################################################################################
library(bayeshscore)
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
# Case 3: true model = AR(1)
#############################################################################################
true_model = 1
true_theta = c(0.6,0.8)
observations1 = simulateData(model(true_model),true_theta,nobservations)
# observations in a matrix of dimensions dimy x nobservations
#--------------------------------------------------------------------------------------------
results_all1 = data.frame()
#--------------------------------------------------------------------------------------------
### Compute logevidence and hscore for each model
for (m in 1:nb_models){
  results = foreach(i=1:repl,.packages=c('bayeshscore'),.verbose = TRUE) %dorng% {
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
# Case 4: true model = AR(2)
#############################################################################################
true_model = 2
true_theta = c(0.25,0.5,0.75^2)
observations2 = simulateData(model(true_model),true_theta,nobservations)
# observations in a matrix of dimensions dimy x nobservations
#--------------------------------------------------------------------------------------------
results_all2 = data.frame()
#--------------------------------------------------------------------------------------------
### Compute logevidence and hscore for each model
for (m in 1:nb_models){
  results = foreach(i=1:repl,.packages=c('bayeshscore'),.verbose = TRUE) %dorng% {
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
  ylab("log-Bayes factor  [1 vs 2]") + facet_grid(. ~ type) + xlab("Number of observations") +
  # guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave("example_AR1_AR2_logBF_1_vs_2.png",dpi = 300,width = 10,height = 5)


# Hyvarinen factor
labels.df = data.frame(x = rep(375,2), y = c(-220,30),
                       text = c("Case 4","Case 3"),
                       type = factor(c("Case 4","Case 3")))
colors = c(wes_palette("Darjeeling")[1],"mediumblue")
axis_ticktextsize = 10
axis_titlesize = 12
ggplot() +
  geom_label(data = labels.df, aes(x,y,label = text,color=type), color = colors, fontface = "bold") +
  theme(legend.position="none") +
  scale_color_manual(values = colors[2:1]) +
  geom_line(data = h_factors, aes(time, hfactor, color = case, group = interaction(case,repl)),alpha=0.6) +
  geom_hline(yintercept = 0,linetype=2) +
  xlab("Number of observations") +
  ylab("H-factor  [1 vs. 2]") +
  ylim(c(min(h_factors$hfactor),30+max(h_factors$hfactor))) +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(size = axis_ticktextsize),
        axis.text.y = element_text(size = axis_ticktextsize),
        axis.title.x = element_text(size = axis_titlesize, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = axis_titlesize, angle = 90, margin = margin(0,10,0,0)),
        strip.text.x = element_text(size = axis_titlesize, colour = "black"),
        strip.background = element_rect(fill="gray88"),
        panel.background = element_rect(fill="gray95",linetype = "solid", colour="white"),
        legend.position = "none")
# ggsave("example_AR1_AR2_Hyvarinen_factor_1_vs_2_10_by_5.png",dpi = 300,width = 10,height = 5)
# ggsave("example_AR1_AR2_Hyvarinen_factor_1_vs_2_5_by_5.png",dpi = 300,width = 5,height = 5)


# log Bayes factor
labels.df = data.frame(x = rep(375,2), y = c(-70,12),
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
  ylab("log-Bayes factor  [1 vs. 2]") +
  ylim(c(min(h_factors$hfactor),30+max(h_factors$hfactor))) +
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
# ggsave("example_AR1_AR2_logBF_1_vs_2_10_by_5.png",dpi = 300,width = 10,height = 5)
# ggsave("example_AR1_AR2_logBF_1_vs_2_5_by_5.png",dpi = 300,width = 5,height = 5)

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
crit.df = data.frame()
crit.df = rbind(crit.df,data.frame(logbayesfactors[-4], value = logbayesfactors$logbayesfactor, crit = factor("LBF")))
crit.df = rbind(crit.df,data.frame(h_factors[-4], value = h_factors$hfactor, crit = factor("HF")))
relabelling = function(x){
  if(x==levels(crit.df$type[1])[1]){return (factor(3))}
  if(x==levels(crit.df$type[1])[2]){return (factor(4))}
}
crit.df$type = sapply(crit.df$type, function(x)relabelling(x))
#--------------------------------------------------------------------------------------------
case_label <- list(
  '3'=expression(paste("Case 3: ",M[1]," nested in ", M[2],", both well-specified",sep="")),
  '4'=expression(paste("Case 4: ",M[1]," nested in ", M[2], ", only ", M[2]," well-specified",sep=""))
)
case_labeller <- function(variable,value){
  return(case_label[value])
}
labels.df = data.frame(x = c(900,900,rep(875,2)),
                       y = c(5.5,1.5,-280,-65),
                       text = rep(c("HF 1 vs. 2","log-BF 1 vs. 2"),4),
                       crit = factor(rep(c("HF","LBF"),4)),
                       type = factor(rep(3:4,each=2)))
hline.df = data.frame(y = c(-5,10), type = factor(3))
axis_titlesize = 18
axis_ticktextsize = 15
ggplot() +
  geom_label(data = labels.df, aes(x,y,label = text,color=crit), fontface = "bold",size=5) +
  theme(legend.position="none") +
  scale_color_manual(values = colors[2:1]) +
  geom_line(data = crit.df, aes(time, value, color = crit, group = interaction(crit,repl)),alpha=0.6) +
  geom_hline(yintercept = 0,linetype=2) +
  geom_hline(data = hline.df, aes(yintercept = y),linetype=2, alpha = 0) +
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
        strip.background = element_rect(fill="gray88"),
        panel.background = element_rect(fill="gray95",linetype = "solid", colour="white"),
        legend.position = "none")
# ggsave("example_AR1_AR2_15_by_6.png",dpi = 300,width = 15,height = 6)
# ggsave("example_AR1_AR2_15_by_6.pdf",dpi = 300,width = 15,height = 6)
