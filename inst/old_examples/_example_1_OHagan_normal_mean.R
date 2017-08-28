##################################################################################################
# This implements example 3.1. in O'Hagan (1995).
# The model is N(theta,1) with conjugate prior on theta, i.e. theta follows N(0,sigma2prior).
# We compute the logevidence and the prequential Hyvarinen score for increasing vagueness
# (i.e. increasing values of sigma2prior).
##################################################################################################
rm(list = ls())
library(HyvarinenSSM)
library(ggplot2)
library(gridExtra)
library(doParallel)
library(foreach)
library(wesanderson)
set.seed(29)

# Define model and data
nobservations = 30
Y = rnorm(nobservations,0,1)
observations = matrix(Y, nrow = 1)# observations in a matrix of dimensions dimy x nobservations

#--------------------------------------------------------------------------------------------
# set algorithmic parameters
algorithmic_parameters = list()
algorithmic_parameters$Ntheta = 2^10
algorithmic_parameters$verbose = TRUE
algorithmic_parameters$ess_threshold = 0.5
# The remaining algorithmic parameters are set to their default values via the functions in util_default.R
#--------------------------------------------------------------------------------------------
repl = 50 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
#--------------------------------------------------------------------------------------------
sigma2prior_all = c(100,1000,10000,100000)
results_all = data.frame()
post_all = data.frame()
#--------------------------------------------------------------------------------------------
### Compute logevidence and hscore for increasing vagueness
for (s in 1:length(sigma2prior_all)){
  sigma2prior = sigma2prior_all[s]
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations, get_model_iid_gaussian_unknown_mean(0,sigma2prior), algorithmic_parameters)
  }
  for (r in 1:repl){
    results_all = rbind(results_all,data.frame(logevidence = results[[r]]$logevidence,
                                               hscore = results[[r]]$hscore,
                                               time = 1:nobservations,
                                               sigma2prior = sigma2prior,
                                               repl = r))
    post_all = rbind(post_all,data.frame(theta = c(results[[r]]$thetas),
                                                   W = c(results[[r]]$normw),
                                                   sigma2prior = sigma2prior,
                                                   repl = r))
  }
}
colors = wes_palette("GrandBudapest")[c(1,4,2,3)]
#--------------------------------------------------------------------------------------------
# Checking sample from the posterior distribution (marginal histogram)
ggplot(post_all, aes(color=factor(format(sigma2prior, scientific = FALSE)))) +
  geom_density(aes(theta,weight=W,group=interaction(sigma2prior,repl))) +
  # scale_color_discrete(expression(paste(" ",sigma[0]^2))) +
  scale_color_manual(expression(bold(paste(" ",sigma[0]^2))),values = colors) +
  xlab(expression(theta)) + guides(colour = guide_legend(override.aes = list(size=2)))
#--------------------------------------------------------------------------------------------
# Check the log-evidence
ggplot(results_all, aes(color=factor(format(sigma2prior, scientific = FALSE)), shape = factor(format(sigma2prior, scientific = FALSE)))) +
  geom_line(aes(time, -logevidence, group=interaction(sigma2prior,repl))) +
  # geom_point(aes(time, -logevidence)) +
  # scale_colour_discrete(expression(paste(" ",sigma[0]^2))) +
  scale_color_manual(expression(bold(paste(" ",sigma[0]^2))),values = colors) +
  # scale_shape_manual("", values=15:18) +
  # scale_colour_discrete("") +
  ylab("- log evidence") + guides(colour = guide_legend(override.aes = list(size = 2)))
#--------------------------------------------------------------------------------------------
# Check the h-score
ggplot(results_all) +
  geom_line(aes(time, hscore, color = factor(format(sigma2prior, scientific = FALSE)),group=interaction(sigma2prior,repl))) +
  # scale_colour_discrete(expression(paste(" ",sigma[0]^2))) +
  scale_color_manual(expression(bold(paste(" ",sigma[0]^2))),values = colors) +
  ylab("Hyvarinen score") + guides(colour = guide_legend(override.aes = list(size=2)))
#--------------------------------------------------------------------------------------------
# Check the log-evidence
ggplot(subset(results_all,time == nobservations)) +
  scale_color_manual(expression(bold(paste(" ",sigma[0]^2))),values = colors) +
  geom_boxplot(aes(factor(sigma2prior),-logevidence, color = factor(sigma2prior)),outlier.shape = NA,size=1) +
  ylab("- log evidence") + xlab(expression(paste(" ",sigma[0]^2))) + theme(legend.position="none")
#--------------------------------------------------------------------------------------------
# Check the h-score
ggplot(subset(results_all,time == nobservations)) +
  geom_boxplot(aes(factor(sigma2prior),hscore, color = factor(sigma2prior)),outlier.shape = NA,size=1) +
  scale_color_manual(expression(bold(paste(" ",sigma[0]^2))),values = colors) +
  ylab("Hyvarinen score") + xlab(expression(paste(" ",sigma[0]^2))) + theme(legend.position="none")
# #--------------------------------------------------------------------------------------------
# # Check the partial log-evidence with training sample size m = 1
# m = 1
# partial_bayes = data.frame()
# for (i in 1: length(sigma2prior_all)){
#   for (r in 1:repl){
#     logevidence = subset(results_all,sigma2prior== sigma2prior_all[i]&repl==r)$logevidence
#     partiallogevidence = logevidence[(m+1):nobservations] - logevidence[1:m]
#     partial_bayes = rbind(partial_bayes,data.frame(time = (m+1):nobservations, sigma2prior = factor(format(sigma2prior_all[i], scientific = FALSE)), repl = r,
#                                                    partiallogevidence = partiallogevidence))
#   }
# }
# ggplot(partial_bayes) +
#   geom_line(aes(time, -partiallogevidence, color = sigma2prior,group=interaction(sigma2prior,repl))) +
#   scale_colour_discrete(expression(paste(" ",sigma[0]^2)))

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Generate plots for paper
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
criteria.df = data.frame()
for (s in 1:length(sigma2prior_all)){
  sigma2 = sigma2prior_all[s]
  for (r in 1:repl){
    result = subset(results_all,sigma2prior==sigma2&repl==r)
    criteria.df = rbind(criteria.df,data.frame(time = 1:nobservations,
                                               repl = r,
                                               sigma2prior = sigma2,
                                               value = -result$logevidence,
                                               type = factor("- log evidence")))
    criteria.df = rbind(criteria.df,data.frame(time = 1:nobservations,
                                               repl = r,
                                               sigma2prior = sigma2,
                                               value = result$hscore,
                                               type = factor("Hyvärinen score")))
  }
}

# plot logevidence and hscore
ggplot(criteria.df, aes(color=factor(format(sigma2prior, scientific = FALSE)), group = interaction(type,repl,sigma2prior))) +
  geom_line(aes(time,value)) +
  # scale_color_manual(expression(bold(paste(" ",sigma[0]^2))),values = wes_palette("Zissou")[c(1,2,4,5)]) +
  ylab("") + xlab("Number of observations") + facet_grid(type ~ ., scales="free") +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0))) +
  scale_color_manual(expression(bold(paste(" ",sigma[0]^2))),values = colors)
# ggsave("example_1_OHagan_normal_mean_1.png",dpi = 300,width = 10,height = 5)



# boxplot of final logevidence and hscore
scaleFUNy <- function(x) sprintf("%.1f", x)
equal_breaks <- function(n, s, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    seq(min(x)+d, max(x)-d, length=n)
  }
}
ggplot(subset(criteria.df,time == nobservations)) +
  scale_color_manual(expression(bold(paste(" ",sigma[0]^2))),values = colors) +
  geom_boxplot(aes(factor(format(sigma2prior, scientific = FALSE)),value, color = factor(format(sigma2prior, scientific = FALSE))),outlier.shape = NA,size=1) +
  xlab(expression(paste(" ",sigma[0]^2))) +
  facet_grid(type ~ ., scales="free") + ylab("") +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0))) +
  scale_y_continuous(labels=scaleFUNy,breaks=equal_breaks(n=5, s=0.1))
# ggsave("example_1_OHagan_normal_mean_2.png",dpi = 300,width = 10,height = 5)
#--------------------------------------------------------------------------------------------

ggplot(subset(criteria.df,time == nobservations)) +
  scale_color_manual(expression(bold(paste(" ",sigma[0]^2))),values = colors) +
  geom_boxplot(aes(factor(format(sigma2prior, scientific = FALSE)),value, color = factor(format(sigma2prior, scientific = FALSE))),outlier.shape = NA,size=1) +
  xlab(expression(paste(" ",sigma[0]^2))) +
  facet_grid(type ~ ., scales="free") + ylab("") +
  theme(strip.text.y = element_text(size = 16, colour = "black", face="plain")) +
  theme(legend.position="none") +
  theme(legend.text=element_text(size=20), legend.title=element_text(size=20)) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin = margin(0,10,0,0))) +
  scale_y_continuous(labels=scaleFUNy,breaks=equal_breaks(n=5, s=0.1))
# ggsave("poster_example_1_OHagan_normal_mean_2.png",dpi = 300,width = 10,height = 5)
