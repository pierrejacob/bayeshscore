##################################################################################################
# Example - iid Normal (example 3.2. in O'Hagan, 1995)
##################################################################################################
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(wesanderson)
set.seed(19)
#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$store_theta = TRUE
#--------------------------------------------------------------------------------------------
# set hyperparameters
muprior = 0
sigma2prior = 100
nu0 = 0.1
s02 = 1
# define models
model = function(i){
  if(i==1){return(get_model_iid_gaussian_unknown_mean(muprior,sigma2prior))} #iid N(theta1, 1)
  if(i==2){return(get_model_iid_gaussian_unknown_variance(nu0,s02))} #iid N(0, theta2)
}
#--------------------------------------------------------------------------------------------
nobservations = 1000
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------
#############################################################################################
# Case 1: DGP = N(1,1), Model 1 is well-specified
# Case 2: DGP = N(0,5), Model 2 is well-specified
# Case 3: DGP = N(2,3), both model 1 and 2 are misspecified
# Case 4: DGP = N(0,1), both model 1 and 2 are well-specified
#############################################################################################
DGP_mu = c(1,0,2,0)
DGP_sigma2 = c(1,5,3,1)
#--------------------------------------------------------------------------------------------
results_all = list()
for (i in 1:length(DGP_mu)){
  results_all_i = data.frame()
  # Generate observations
  observations = matrix(rnorm(nobservations,DGP_mu[i],sqrt(DGP_sigma2[i])), nrow = 1)
  # Compute logevidence and hscore
  for (m in 1:2){
    results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
      hscore(observations, model(m), algorithmic_parameters)
    }
    for (r in 1:repl){
      results_all_i = rbind(results_all_i,data.frame(logevidence = results[[r]]$logevidence,
                                                     hscore = results[[r]]$hscore,
                                                     time = 1:nobservations,
                                                     model = factor(m),
                                                     repl = factor(r)))
    }
  }
  results_all[[i]] = results_all_i
}
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Generate plots for paper
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Compute the Hyvarinen factor
h_factors = data.frame()
plot_hfactor = list()
for (i in 1:4){
  for (r in 1:repl) {
    h_factor = subset(results_all[[i]],model==2&repl==r)$hscore-subset(results_all[[i]],model==1&repl==r)$hscore
    h_factors = rbind(h_factors,data.frame(time = 1:nobservations,
                                           repl = r,
                                           hfactor = h_factor,
                                           case = factor(i),
                                           type = factor(paste("Case",toString(i))),
                                           sim = 1))
  }
}
# Checking H-factor
# top-left, top-right, bottom-left, bottom-right = case 1, 2, 3, 4.
# Positive = choose model 1 // Negative == choose model 2.
case_label <- list(
  'Case 1'=expression(paste("Case 1: ",M[1]," is well-specified",sep="")),
  'Case 2'=expression(paste("Case 2: ",M[2]," is well-specified",sep="")),
  'Case 3'=expression(paste("Case 3: both are misspecified",sep="")),
  'Case 4'=expression(paste("Case 4: both are well-specified",sep=""))
)
case_labeller <- function(variable,value){
  return(case_label[value])
}
colors = c("dodgerblue")
ggplot(h_factors, aes(color = factor(sim), group = interaction(case,repl), linetype = factor(sim))) +
  geom_line(aes(time, hfactor),alpha=0.6) +
  # scale_linetype_manual(values=c("dashed","solid")) +
  scale_color_manual(values=colors) +
  geom_hline(yintercept = 0,linetype = 2) +
  xlab("Number of observations") +
  ylab("Hyvärinen factor  [1 vs 2]") +
  # theme_bw() +
  facet_wrap( ~ type, ncol=2, scales="free", labeller = case_labeller) +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave("example_consistency_iidNormal.png",dpi = 300,width = 10,height = 5)


# Compute the Bayes factor
logBFs = data.frame()
plot_logBF = list()
for (i in 1:4){
  for (r in 1:repl) {
    logBF = -subset(results_all[[i]],model==2&repl==r)$logevidence+subset(results_all[[i]],model==1&repl==r)$logevidence
    logBFs = rbind(logBFs,data.frame(time = 1:nobservations,
                                           repl = r,
                                           logBF = logBF,
                                           case = factor(i),
                                           type = factor(paste("Case",toString(i))),
                                           sim = 1))
  }
}
# Checking H-factor
# top-left, top-right, bottom-left, bottom-right = case 1, 2, 3, 4.
# Positive = choose model 1 // Negative == choose model 2.
case_label <- list(
  'Case 1'=expression(paste("Case 1: ",M[1]," is well-specified",sep="")),
  'Case 2'=expression(paste("Case 2: ",M[2]," is well-specified",sep="")),
  'Case 3'=expression(paste("Case 3: both are misspecified",sep="")),
  'Case 4'=expression(paste("Case 4: both are well-specified",sep=""))
)
case_labeller <- function(variable,value){
  return(case_label[value])
}
colors = c("tomato3")
ggplot(logBFs, aes(color = factor(sim), group = interaction(case,repl), linetype = factor(sim))) +
  geom_line(aes(time, logBF),alpha=0.6) +
  # scale_linetype_manual(values=c("dashed","solid")) +
  scale_color_manual(values=colors) +
  geom_hline(yintercept = 0,linetype = 2) +
  xlab("Number of observations") +
  ylab("log Bayes factor  [1 vs 2]") +
  # theme_bw() +
  facet_wrap( ~ type, ncol=2, scales="free", labeller = case_labeller) +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave("example_consistency_iidNormal_logBF.png",dpi = 300,width = 10,height = 5)


# ################################################################################################
# ################################################################################################
# # Plot for poster
# ################################################################################################
# ################################################################################################
# ggplot(h_factors, aes(color = factor(sim), group = interaction(case,repl), linetype = factor(sim))) +
#   geom_line(aes(time, hfactor)) +
#   # scale_linetype_manual(values=c("dashed","solid")) +
#   scale_color_manual(values=colors) +
#   geom_hline(yintercept = 0,alpha=0.3) +
#   ylab("Hyvärinen factor  [1 vs 2]") +
#   facet_wrap( ~ type, ncol=2, scales="free", labeller = case_labeller) +
#   xlab("Number of observations") +
#   theme(strip.text.x = element_text(size = 16, colour = "black", face="plain")) +
#   theme(legend.position="none") +
#   theme(legend.text=element_text(size=20), legend.title=element_text(size=20)) +
#   theme(axis.text.x = element_text(size = 15),
#         axis.text.y = element_text(size = 15),
#         axis.title.x = element_text(size = 20, margin = margin(10,0,0,0)),
#         axis.title.y = element_text(size = 20, margin = margin(0,10,0,0)))
# # ggsave("poster_example_2_OHagan_Hyvarinen_factor_1_versus_2.png",dpi = 300,width = 10,height = 5)

