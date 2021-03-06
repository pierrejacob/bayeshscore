##################################################################################################
# Example - iid Normal (example 3.2. in O'Hagan, 1995)
##################################################################################################
library(bayeshscore)
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
algorithmic_parameters$store_theta = FALSE
#--------------------------------------------------------------------------------------------
# set hyperparameters
muprior = 0
sigma2prior = 10
nu0 = 0.1
s02 = 1
# define models
model = function(i){
  if(i==1){return(get_model_iid_gaussian_unknown_mean(muprior,sigma2prior))} #iid N(theta1, 1)
  if(i==2){return(get_model_iid_gaussian(muprior,sigma2prior,nu0,s02))} #iid N(0, theta2)
}
models_to_run = c(1,2)
#--------------------------------------------------------------------------------------------
nobservations = 1000
M = 100
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
#############################################################################################
# Case 2: DGP = N(0,5), Model 2 is well-specified
# Case 4: DGP = N(0,1), both model 1 and 2 are well-specified
#############################################################################################
DGP_mu = c(0,0)
DGP_sigma2 = c(5,1)
#--------------------------------------------------------------------------------------------
observations = list()
for (i in 1:length(DGP_mu)){
  # Generate observations
  observations[[i]] = matrix(rnorm(nobservations,DGP_mu[i],sqrt(DGP_sigma2[i])), nrow = 1)
}
# Run SMC
results_all = list()
post_all = list()
schedule = expand.grid(1:repl,1:2,1:length(DGP_mu))
results = foreach(s=1:nrow(schedule),.packages=c('bayeshscore'),.verbose = TRUE) %dorng% {
  dgp = schedule[s,3]
  m = schedule[s,2]
  r = schedule[s,1]
  sink(logfilename, append = TRUE) # Monitor progress in parallel via log file
  average_hscore = 0
  average_logevidence = 0
  for (permut in 1:M){
    observations_permuted = observations[[dgp]][,sample(nobservations),drop=FALSE]
    hscore_perm = hscore(observations_permuted, model(m), algorithmic_parameters)
    average_hscore = average_hscore + hscore_perm$hscore
    average_logevidence = average_logevidence + hscore_perm$logevidence
  }
  c(hscore_perm, list(avrg_hscore = average_hscore/M,avrg_logevidence = average_logevidence/M))
}
for (dgp in 1:length(DGP_mu)){
  for (m in models_to_run) {
    for (r in 1:repl) {
      s = (1:nrow(schedule))[(schedule[,1]==r) & (schedule[,2]==m) & (schedule[,3]==dgp)]
      results_all = rbind(results_all,data.frame(logevidence = results[[s]]$avrg_logevidence,
                                                 hscore = results[[s]]$avrg_hscore,
                                                 time = 1:nobservations,
                                                 model = factor(m),
                                                 repl = factor(r),
                                                 case = factor(dgp)))
      # post_all = rbind(post_all,data.frame(theta = t(results[[s]]$thetas),
      #                                      W = results[[s]]$normw,
      #                                      repl = factor(r),
      #                                      model = factor(m),
      #                                      case = factor(dgp)))
    }
  }
}
#--------------------------------------------------------------------------------------------
# close log file
sink()
#--------------------------------------------------------------------------------------------
# Checking posterior distributions (marginal)
#--------------------------------------------------------------------------------------------
#Compute exact posterior
# sigma2_post = rep(NA, length(DGP_mu))
# mu_post = rep(NA, length(DGP_mu))
# nu_post = rep(NA, length(DGP_mu))
# s2_post = rep(NA, length(DGP_mu))
# plot_post = list()
# for (i in 1:length(DGP_mu)) {
#   sigma2_post[i] = 1/(nobservations + 1/sigma2prior)
#   mu_post[i] = (sum(observations[[i]]) + (1/sigma2prior)*muprior)*sigma2_post[i]
#   nu_post[i] = nu0 + nobservations
#   s2_post[i] = (nu0*s02 + sum(observations[[i]]^2))/nu_post[i]
#   # Plot posterior samples vs exact posterior
#   local({i = i;
#   # model 1
#   plot_post[[2*i-1]] <<- ggplot(subset(post_all,case==i & model==1)) +
#     geom_density(aes(theta, weight = W, fill = factor(i), group = interaction(repl,model)), alpha = 0.3) +
#     stat_function(fun = function(y)dnorm(y,mu_post[i],sqrt(sigma2_post[i]),FALSE),colour="blue",size=1.5,linetype=1) +
#     theme(legend.position="none") + xlab(expression(theta[1])) + ylab("");
#   # model 2
#   plot_post[[2*i]] <<- ggplot(subset(post_all,case==i & model==2)) +
#     geom_density(aes(theta, weight = W, fill = factor(i), group = interaction(repl,model)), alpha = 0.3) +
#     stat_function(fun = function(y)dinvchisq(y,nu_post[i],s2_post[i],FALSE),colour="blue",size=1.5,linetype=1) +
#     xlab(expression(theta[2])) + ylab("") + scale_fill_discrete(name = "case");
#   })
# }
# # Checking posterior distributions (marginal)
# # One row for each case: theta1 for model 1, theta2 for model 2
# grid.arrange(plot_post[[1]],plot_post[[2]],
#              plot_post[[3]],plot_post[[4]],
#              plot_post[[5]],plot_post[[6]],
#              plot_post[[7]],plot_post[[8]],
#              ncol = 4, widths = c(1,1.4,1,1.4))
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Generate plots for paper
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Compute the Hyvarinen factor and log-Bayes factor
criteria.df = data.frame()
for (i in 1:length(DGP_mu)){
  for (r in 1:repl) {
    m1 = subset(results_all,case==i&model==1&repl==r)
    m2 = subset(results_all,case==i&model==2&repl==r)
    h_factor = m2$hscore - m1$hscore
    logbayes_factor = -m2$logevidence + m1$logevidence
    criteria.df = rbind(criteria.df,data.frame(time = 1:nobservations,
                                               repl = factor(r),
                                               value = h_factor,
                                               case = factor(i),
                                               crit = factor("HF"),
                                               sim = factor(1)))
    criteria.df = rbind(criteria.df,data.frame(time = 1:nobservations,
                                               repl = factor(r),
                                               value = logbayes_factor,
                                               case = factor(i),
                                               crit = factor("LBF"),
                                               sim = factor(1)))
  }
}
# Compute theoretical slopes
getDivergence = function(mu_star,sigma2_star){
  DH1 = (sigma2_star-1)^2/sigma2_star
  DH2 = (mu_star/(mu_star^2+sigma2_star))^2*(1+mu_star^2/sigma2_star)
  KL1 = (sigma2_star-1-log(sigma2_star))/2
  KL2 = log(1+mu_star^2/sigma2_star)/2
  return (list(DH1 = DH1, DH2 = DH2, KL1 = KL1, KL2 = KL2))
}
slopes.df = data.frame()
for (i in 1:length(DGP_mu)){
  divs = getDivergence(DGP_mu[i],DGP_sigma2[i])
  slopes.df = rbind(slopes.df,data.frame(time = 1:nobservations,
                                         slope = (1:nobservations)*(divs$DH2 - divs$DH1),
                                         case = factor(i),
                                         crit = factor("HF"),
                                         sim = factor(0)))
  slopes.df = rbind(slopes.df,data.frame(time = 1:nobservations,
                                         slope = (1:nobservations)*(divs$KL2 - divs$KL1),
                                         case = factor(i),
                                         crit = factor("LBF"),
                                         sim = factor(0)))
}
# Plot the H-factor and log-Bayes factor
# top-left, top-right, bottom-left, bottom-right = case 1, 2, 3, 4.
# Positive = choose model 1 // Negative == choose model 2.
case_label <- list(
  '1'=expression(paste("Case 1: ",M[1]," misspecified, ", M[2]," well-specified",sep="")),
  '2'=expression(paste("Case 2: ",M[1]," and ", M[2]," both well-specified",sep=""))
)
case_labeller <- function(variable,value){
  return(case_label[value])
}
labels.df = data.frame(x = c(850,850,900,875),
                       y = c(-2000,-500,10.5,6.9),
                       text = rep(c("HF 1 vs. 2","log-BF 1 vs. 2"),4),
                       crit = factor(rep(c("HF","LBF"),4)),
                       case = factor(rep(1:2,each=2)))
# compute theoretical log-rate
slopes.df$slope[slopes.df$case==2 & slopes.df$crit=="HF"] = 2*log(1:nobservations)
slopes.df$slope[slopes.df$case==2 & slopes.df$crit=="LBF"] = (1/2)*log(1:nobservations)


axis_titlesize = 22
axis_ticktextsize = 15
colors = c("dodgerblue","tomato3")
# plot one point every xstep observations (to speed up the printing and display of the graphs)
x_step = 1
# plot results
ggplot() +
  # scale_linetype_manual(values=c("dashed","solid")) +
  scale_color_manual(values=colors) +
  geom_hline(yintercept = 0,linetype = 2) +
  geom_label(data = labels.df, aes(x,y,label = text,color=crit,fontface="bold"),size=5) +
  geom_line(data = subset(slopes.df,time%in%seq(1,nobservations, x_step)),
            aes(time, slope, group = interaction(crit,case), color = crit), linetype=4, alpha = 1, size = 1.5) +
  geom_line(data = subset(criteria.df,time%in%seq(1,nobservations, x_step)), aes(time, value, color = crit, group = interaction(case,repl,crit)),alpha=0.7) +
  xlab("Number of observations") +
  ylab("") +
  facet_wrap( ~ case, ncol=2, scales="free", labeller = case_labeller) +
  theme(axis.text.x = element_text(size = axis_ticktextsize),
        axis.text.y = element_text(size = axis_ticktextsize),
        axis.title.x = element_text(size = axis_titlesize, margin=margin(20,0,0,0)),
        axis.title.y = element_text(size = axis_titlesize, angle = 90, margin = margin(0,20,0,0)),
        strip.text.x = element_text(size = axis_titlesize, colour = "black"),
        strip.background = element_rect(fill="gray88"),
        panel.background = element_rect(fill="gray95",linetype = "solid", colour="white"),
        legend.position = "none")
# ggsave("example_iidNormal_nested_15_by_9_T=1000_M=100.png",dpi = 300,width = 15,height = 9)
# ggsave("example_iidNormal_nested_15_by_6_T=1000_M=100.pdf",dpi = 300,width = 15,height = 6)
