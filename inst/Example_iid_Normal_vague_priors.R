##################################################################################################
# Example - iid Normal (example 3.2. in O'Hagan, 1995)
##################################################################################################
library(bayeshscore)
library(ggplot2)
library(gridExtra)
library(wesanderson)
library(numDeriv)
set.seed(19)
#--------------------------------------------------------------------------------------------
# set hyperparameters
muprior = 0
logsigmapriors = c(0,150,350)
nu0 = 0.1
s02 = 1
#--------------------------------------------------------------------------------------------
nobservations = 1000
#--------------------------------------------------------------------------------------------
#############################################################################################
# Case 1: DGP = N(1,1), Model 1 is well-specified
#############################################################################################
results = data.frame()
M = 100
####
sigma2_post = array(NA,c(nobservations,length(logsigmapriors),M))
mu_post = array(NA,c(nobservations,length(logsigmapriors),M))
nu_post = array(NA,c(nobservations,length(logsigmapriors),M))
s2_post = array(NA,c(nobservations,length(logsigmapriors),M))
####
pb = txtProgressBar(0,M*length(logsigmapriors),1,style = 3)
count = 0
for (r in 1:M){
  Y = matrix(rnorm(nobservations,1,1), nrow = 1)
  for (j in 1:length(logsigmapriors)){
    sigma2prior = exp(2*logsigmapriors[j])
    logscore1 = rep(NA, nobservations)
    logscore2 = rep(NA, nobservations)
    hscore1 = rep(NA, nobservations)
    hscore2 = rep(NA, nobservations)
    for (t in 1:nobservations) {
      # Compute exact posterior densities for models 1 and 2
      sigma2_post[t,j,r] = 1/(t + 1/sigma2prior)
      mu_post[t,j,r] = (sum(Y[,1:t]) + (1/sigma2prior)*muprior)*sigma2_post[t,j,r]
      nu_post[t,j,r] = nu0 + t
      s2_post[t,j,r] = (nu0*s02 + sum(Y[,1:t]^2))/nu_post[t,j,r]
      # Compute exact predictive log-densities for models 1 and 2
      logd1 = function(y){dnorm(y,ifelse(t==1,muprior,mu_post[t-1,j,r]),sqrt(ifelse(t==1,sigma2prior,1+sigma2_post[t-1,j,r])),TRUE)}
      logd2 = function(y){dtscaled(y,ifelse(t==1,nu0,nu_post[t-1,j,r]),ifelse(t==1,s02,s2_post[t-1,j,r]),TRUE)}
      # Compute exact incremental log-score for models 1 and 2
      logscore1[t] = -logd1(Y[,t])
      logscore2[t] = -logd2(Y[,t])
      # Compute exact incremental h-score for models 1 and 2
      hscore1[t] = 2*hessian(logd1,Y[,t]) + (grad(logd1,Y[,t]))^2
      hscore2[t] = 2*hessian(logd2,Y[,t]) + (grad(logd2,Y[,t]))^2
    }
    # Compute exact log-score for models 1 and 2
    logscore1 = cumsum(logscore1)
    logscore2 = cumsum(logscore2)
    # Compute exact h-score for models 1 and 2
    hscore1 = cumsum(hscore1)
    hscore2 = cumsum(hscore2)
    # Format the results into a dataframe
    results = rbind(results,data.frame(value = c(logscore2-logscore1, hscore2-hscore1),
                                       type = factor(rep(c("LBF","HF"),each = nobservations)),
                                       time = rep(1:nobservations,2),
                                       repl = factor(r,levels = 1:M),
                                       case = factor(j,levels = 1:length(logsigmapriors))))
    count = count + 1
    setTxtProgressBar(pb,count)
  }
}
axis_titlesize = 22
axis_ticktextsize = 15
colors = c("dodgerblue","tomato3")
case_label <- list(
  '1'=expression(paste("log(",sigma[0],") = 0",sep="")),
  '2'=expression(paste("log(",sigma[0],") = 150",sep="")),
  '3'=expression(paste("log(",sigma[0],") = 350",sep=""))
)
case_labeller <- function(variable,value){
  return(case_label[value])
}
labels.df = data.frame(x = c(250,750,250,750,250,750),
                       y = c(400,100,400,-100,400,-250),
                       text = rep(c("HF 1 vs. 2","log-BF 1 vs. 2"),3),
                       crit = factor(rep(c("HF","LBF"),3)),
                       case = factor(rep(1:3,each=2)))
# coarsen the curves (i.e. only plot 1 every "time_step" observations) to reduce the figure file size
time_step = 10
ggplot(subset(results, time %in% seq(1,nobservations,time_step))) +
  geom_hline(yintercept = 0,linetype = 2) +
  scale_color_manual(values=colors) +
  geom_line(aes(time, value, group = interaction(repl,type), color = type), alpha = 0.3) +
  geom_label(data = labels.df, aes(x,y,label = text,color=crit,fontface="bold"),size=6.5) +
  facet_wrap( ~ case, ncol=3, labeller = case_labeller) +
  xlab("Number of observations") +  ylab("") +
  theme(axis.text.x = element_text(size = axis_ticktextsize),
        axis.text.y = element_text(size = axis_ticktextsize),
        axis.title.x = element_text(size = axis_titlesize, margin=margin(20,0,0,0)),
        axis.title.y = element_text(size = axis_titlesize, angle = 90, margin = margin(0,20,0,0)),
        strip.text.x = element_text(size = axis_titlesize, colour = "black"),
        strip.background = element_rect(fill="gray88"),
        panel.background = element_rect(fill="gray95",linetype = "solid", colour="white"),
        legend.position = "none")
# ggsave("example_iidNormal_vague_15_by_7.5.png",dpi = 300,width = 15,height = 7.5)
# ggsave("example_iidNormal_vague_15_by_7.5.pdf",dpi = 300,width = 15,height = 7.5)

t = nobservations
posteriors = data.frame()
grid_nb_points = 100
mu_grid = seq(min(mu_post[t,,]-5*sqrt(max(sigma2_post[t,,]))), max(mu_post[t,,]+5*sqrt(max(sigma2_post[t,,]))),
              length.out = grid_nb_points)
sigma2_grid = seq(1.4, 2.7, length.out = grid_nb_points)
for (r in 1:M){
  for (j in 1:length(logsigmapriors)){
    posteriors = rbind(posteriors, data.frame(dens = dnorm(mu_grid,mu_post[t,j,r],sqrt(sigma2_post[t,j,r])),
                                              time = t,
                                              grid = mu_grid,
                                              repl = factor(r,levels = 1:M),
                                              case = factor(logsigmapriors[j], levels = logsigmapriors),
                                              param = factor("mu", levels = c("mu","sigma2"))))
    posteriors = rbind(posteriors, data.frame(dens = dinvchisq(sigma2_grid,nu_post[t,j,r],s2_post[t,j,r],FALSE),
                                              time = t,
                                              grid = sigma2_grid,
                                              repl = factor(r,levels = 1:M),
                                              case = factor(logsigmapriors[j], levels = logsigmapriors),
                                              param = factor("sigma2", levels = c("mu","sigma2"))))
  }
}

ggplot(subset(posteriors,param=="mu")) +
  geom_line(aes(grid, dens, group = repl), alpha = 0.3, color = "forestgreen") +
  xlab("") +  ylab("") +
  facet_grid(. ~ case, scales = "free",
             labeller = label_bquote(cols = log(sigma[0]) == .((logsigmapriors)[case]),
                                     rows = paste("Posterior density of ",mu,sep=""))) +
  theme(axis.text.x = element_text(size = axis_ticktextsize),
        axis.text.y = element_text(size = axis_ticktextsize),
        axis.title.x = element_text(size = axis_titlesize, margin=margin(20,0,0,0)),
        axis.title.y = element_text(size = axis_titlesize, angle = 90, margin = margin(0,20,0,0)),
        strip.text.x = element_text(size = axis_titlesize, colour = "black"),
        strip.text.y = element_text(size = axis_titlesize, colour = "black"),
        strip.background = element_rect(fill="gray88"),
        panel.background = element_rect(fill="gray95",linetype = "solid", colour="white"),
        legend.position = "none")
# ggsave("example_post_mu_vague_15_by_7.5.png",dpi = 300,width = 15,height = 7.5)
# ggsave("example_post_mu_vague_15_by_7.5.pdf",dpi = 300,width = 15,height = 7.5)
