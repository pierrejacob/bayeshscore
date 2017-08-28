##################################################################################################
# This implements example 3.2. in O'Hagan (1995).
##################################################################################################
rm(list = ls())
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
nobservations = 100
#--------------------------------------------------------------------------------------------
repl = 5 #number of replications
registerDoParallel(cores=5) #number of workers in parallel
##################################################################################################
# Case 1: DGP = N(1,1), Model 1 is correct
##################################################################################################
observations1 = matrix(rnorm(nobservations,1,sqrt(1)), nrow = 1)# observations (dimy by nobservations matrix)
#-----------------------------------
results1_all = data.frame()
post1_all = data.frame()
#-----------------------------------
# Compute logevidence and hscore
for (m in 1:2){
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations1, model(m), algorithmic_parameters)
  }
  for (r in 1:repl){
    results1_all = rbind(results1_all,data.frame(logevidence = results[[r]]$logevidence,
                                                 hscore = results[[r]]$hscore,
                                                 time = 1:nobservations,
                                                 model = factor(m),
                                                 repl = factor(r)))
    post1_all = rbind(post1_all,data.frame(theta = c(results[[r]]$thetas),
                                           W = results[[r]]$normw,
                                           model = factor(m),
                                           repl = factor(r)))
  }
}
##################################################################################################
# Case 2: DGP = N(0,5), Model 2 is correct
##################################################################################################
observations2 = matrix(rnorm(nobservations,0,5), nrow = 1)# observations (dimy by nobservations matrix)
#-----------------------------------
results2_all = data.frame()
post2_all = data.frame()
#-----------------------------------
# Compute logevidence and hscore
for (m in 1:2){
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations2, model(m), algorithmic_parameters)
  }
  for (r in 1:repl){
    results2_all = rbind(results2_all,data.frame(logevidence = results[[r]]$logevidence,
                                                 hscore = results[[r]]$hscore,
                                                 time = 1:nobservations,
                                                 model = factor(m),
                                                 repl = factor(r)))
    post2_all = rbind(post2_all,data.frame(theta = c(results[[r]]$thetas),
                                           W = results[[r]]$normw,
                                           model = factor(m),
                                           repl = factor(r)))
  }
}
##################################################################################################
# Case 3: DGP = N(2,3), both model 1 and 2 are misspecified
##################################################################################################
observations3 = matrix(rnorm(nobservations,2,sqrt(3)), nrow = 1)# observations (dimy by nobservations matrix)
#-----------------------------------
results3_all = data.frame()
post3_all = data.frame()
#-----------------------------------
# Compute logevidence and hscore
for (m in 1:2){
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations3, model(m), algorithmic_parameters)
  }
  for (r in 1:repl){
    results3_all = rbind(results3_all,data.frame(logevidence = results[[r]]$logevidence,
                                                 hscore = results[[r]]$hscore,
                                                 time = 1:nobservations,
                                                 model = factor(m),
                                                 repl = factor(r)))
    post3_all = rbind(post3_all,data.frame(theta = c(results[[r]]$thetas),
                                           W = results[[r]]$normw,
                                           model = factor(m),
                                           repl = factor(r)))
  }
}
##################################################################################################
# Case 4: DGP = N(0,1), both model 1 and 2 are correct
##################################################################################################
observations4 = matrix(rnorm(nobservations,0,1), nrow = 1)# observations (dimy by nobservations matrix)
#-----------------------------------
results4_all = data.frame()
post4_all = data.frame()
#-----------------------------------
# Compute logevidence and hscore
for (m in 1:2){
  results = foreach(i=1:repl,.packages=c('HyvarinenSSM'),.verbose = TRUE) %dorng% {
    hscore(observations4, model(m), algorithmic_parameters)
  }
  for (r in 1:repl){
    results4_all = rbind(results4_all,data.frame(logevidence = results[[r]]$logevidence,
                                                 hscore = results[[r]]$hscore,
                                                 time = 1:nobservations,
                                                 model = factor(m),
                                                 repl = factor(r)))
    post4_all = rbind(post4_all,data.frame(theta = c(results[[r]]$thetas),
                                           W = results[[r]]$normw,
                                           model = factor(m),
                                           repl = factor(r)))
  }
}
# ##################################################################################################
# ##################################################################################################
# ##################################################################################################
# ##################################################################################################
# ##################################################################################################
# #--------------------------------------------------------------------------------------------
# # Checking posterior distributions (marginal)
# #--------------------------------------------------------------------------------------------
post_all = list(post1_all, post2_all, post3_all, post4_all)
# observations = list(observations1, observations2, observations3, observations4)
# #Compute exact posterior
# sigma2_post = rep(NA, 4)
# mu_post = rep(NA, 4)
# nu_post = rep(NA, 4)
# s2_post = rep(NA, 4)
# plot_post = list()
# for (i in 1:4) {
#   sigma2_post[i] = 1/(nobservations + 1/sigma2prior)
#   mu_post[i] = (sum(observations[[i]]) + (1/sigma2prior)*muprior)*sigma2_post[i]
#   nu_post[i] = nu0 + nobservations
#   s2_post[i] = (nu0*s02 + sum(observations[[i]]^2))/nu_post[i]
#   # Plot posterior samples vs exact posterior
#   local({i = i;
#   # model 1
#   plot_post[[2*i-1]] <<- ggplot(subset(post_all[[i]], model==1)) +
#     geom_density(aes(theta, weight = W, fill = factor(i), group = interaction(repl,model)), alpha = 0.3) +
#     stat_function(fun = function(y)dnorm(y,mu_post[i],sqrt(sigma2_post[i]),FALSE),colour="blue",size=1.5,linetype=1) +
#     theme(legend.position="none") + xlab(expression(theta[1])) + ylab("");
#   # model 2
#   plot_post[[2*i]] <<- ggplot(subset(post_all[[i]], model==2)) +
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
# #--------------------------------------------------------------------------------------------
# # Sanity check (exact computation)
# #--------------------------------------------------------------------------------------------
results_all = list(results1_all, results2_all, results3_all, results4_all)
# M = 1000
# DGP_mu = c(1,0,2,0)
# DGP_sigma2 = c(1,5,3,1)
# plot_exact = list()
# for (i in 1:4){
#   expected_H_model1_exact = rep(0,4)
#   expected_H_model2_exact = rep(0,4)
#   sigma2_post = rep(NA, nobservations)
#   mu_post = rep(NA, nobservations)
#   nu_post = rep(NA, nobservations)
#   s2_post = rep(NA, nobservations)
#   hexact1 = rep(NA, nobservations)
#   hexact2 = rep(NA, nobservations)
#   for (t in 1:nobservations) {
#     sigma2_post[t] = 1/(t + 1/sigma2prior)
#     mu_post[t] = (sum(observations[[i]][,1:t]) + (1/sigma2prior)*muprior)*sigma2_post[t]
#     nu_post[t] = nu0 + t
#     s2_post[t] = (nu0*s02 + sum(observations[[i]][,1:t]^2))/nu_post[t]
#     theta2 = rinvchisq(M,nu_post[t],s2_post[t])
#     if (t == 1){
#       hexact1[t] = -2/(sigma2prior+1) + (observations[[i]][,t]^2)/((sigma2prior+1)^2)
#     } else {
#       hexact1[t] = -2/(sigma2_post[t-1]+1) + ((observations[[i]][,t]-mu_post[t-1])^2)/((sigma2_post[t-1]+1)^2)
#     }
#
#     hexact2[t] = mean((1/(theta2^2))*((observations[[i]][,t]-0)^2-2*theta2)) + var(1/theta2)*(observations[[i]][,t])^2
#   }
#   hexact1 = cumsum(hexact1)
#   hexact2 = cumsum(hexact2)
#   expected_H_model1_exact[i] = (DGP_sigma2[i]+(DGP_mu[i]-mean(theta1))^2-2)
#   expected_H_model2_exact[i] = (1/mean(theta2)^2)*(DGP_sigma2[i]+DGP_mu[i]^2-2*mean(theta2))
#   local({
#     i = i;
#     hexact1 = hexact1;
#     hexact2 = hexact2;
#     expected_H_model1_exact = expected_H_model1_exact;
#     expected_H_model2_exact = expected_H_model2_exact;
#     results_all = results_all;
#     plot_exact[[i]] <<-ggplot() +
#       geom_line(data=results_all[[i]],aes(time,hscore/time,group=interaction(repl,model),col=model)) +
#       geom_line(aes(1:nobservations,expected_H_model1_exact[i])) +
#       geom_line(aes((1:nobservations),expected_H_model2_exact[i])) +
#       geom_point(aes((1:nobservations),hexact1/(1:nobservations)),col="red") +
#       geom_point(aes((1:nobservations),hexact2/(1:nobservations)),col="blue")
#   })
# }
# do.call(grid.arrange,c(plot_exact, ncol = 2))

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Generate plots for paper
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
colors = c("blue",wes_palette("Royal1")[2])
DGP_mu = c(1,0,2,0)
DGP_sigma2 = c(1,5,3,1)
# Compute the Hyvarinen factor
h_factors = data.frame()
theta1_star = rep(0, 4)
theta2_star = rep(0, 4)
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
    # thetas1 = subset(post_all[[i]],model==1&repl==r)
    # thetas2 = subset(post_all[[i]],model==2&repl==r)
    # theta1_star[i] = theta1_star[i] + sum(thetas1$theta*thetas1$W)
    # theta2_star[i] = theta2_star[i] + sum(thetas2$theta*thetas2$W)
  }
  # # Compute the theoretical (asymptotic) speed of convergence
  # theta1_star[i] = theta1_star[i]/repl
  # theta2_star[i] = theta2_star[i]/repl
  # expected_H_model1 = (DGP_sigma2[i]+(DGP_mu[i]-theta1_star[i])^2-2)
  # expected_H_model2 = (1/theta2_star[i]^2)*(DGP_sigma2[i]+DGP_mu[i]^2-2*theta2_star[i])


  expected_H_model1 = DGP_sigma2[i]-2
  expected_H_model2 = -1/(DGP_sigma2[i] + DGP_mu[i]^2)


  slope = expected_H_model2 - expected_H_model1
  # Plot h-factor (hscore model 2 minus model 1)
  # Thus: positive = in favor of model 1, negative = in favor of model 2

  # h_factors = rbind(h_factors,data.frame(time = 1:nobservations,
  #                                        repl = -1,
  #                                        hfactor = (1:nobservations)*slope,
  #                                        case = factor(i),
  #                                        type = factor(paste("Case",toString(i))),
  #                                        sim = -1))
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
ggplot(h_factors, aes(color = factor(sim), group = interaction(case,repl), linetype = factor(sim))) +
  geom_line(aes(time, hfactor)) +
  # scale_linetype_manual(values=c("dashed","solid")) +
  scale_color_manual(values=colors) +
  geom_hline(yintercept = 0,alpha=0.3) +
  xlab("Number of observations") +
  ylab("Hyvrärinen factor  [1 vs 2]") +
  facet_wrap( ~ type, ncol=2, scales="free", labeller = case_labeller) +
  theme(strip.text.y = element_text(size = 12, colour = "black")) +
  theme(legend.text=element_text(size=12)) +
  theme(legend.title=element_text(size=12)) +
  theme(legend.position="none") +
  theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
  theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
# ggsave("example_2_OHagan_Hyvarinen_factor_1_versus_2.png",dpi = 300,width = 10,height = 5)




################################################################################################
################################################################################################
# Plot for poster
################################################################################################
################################################################################################
ggplot(h_factors, aes(color = factor(sim), group = interaction(case,repl), linetype = factor(sim))) +
  geom_line(aes(time, hfactor)) +
  # scale_linetype_manual(values=c("dashed","solid")) +
  scale_color_manual(values=colors) +
  geom_hline(yintercept = 0,alpha=0.3) +
  ylab("Hyvrärinen factor  [1 vs 2]") +
  facet_wrap( ~ type, ncol=2, scales="free", labeller = case_labeller) +
  xlab("Number of observations") +
  theme(strip.text.x = element_text(size = 16, colour = "black", face="plain")) +
  theme(legend.position="none") +
  theme(legend.text=element_text(size=20), legend.title=element_text(size=20)) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 20, margin = margin(0,10,0,0)))
# ggsave("poster_example_2_OHagan_Hyvarinen_factor_1_versus_2.png",dpi = 300,width = 10,height = 5)

