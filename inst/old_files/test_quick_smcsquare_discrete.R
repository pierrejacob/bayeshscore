rm(list = ls())
library(HyvarinenSSM)
library(gridExtra)
library(doRNG)
library(dplyr)
theme_set(theme_classic())
registerDoMC(cores = 12)
set.seed(17)

# Define model and data
observations <- data_kangaroo[c("y1","y2"),]
# observations <- observations[,1:10]

rangeprior = 10
model2 <- get_model_kangarooExponential(rangeprior)
model3 <- get_model_kangarooRandomwalk(rangeprior)
# Define initial proposal for theta (to avoid sampling from vague prior)
# Define algorithmic parameters for each model
algorithmic_parameters <- list(Ntheta = 2^11, Nx = 2^12,
              resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
              progress = TRUE, store = TRUE)



# library(doMC)
# r3 <-    hscore_discrete(observations, model3, algorithmic_parameters)

R <- 20
r2 <- foreach(r = 1:R) %dorng% {
  hscore_discrete(observations, model2, algorithmic_parameters)
}
r3 <- foreach(r = 1:R) %dorng% {
   hscore_discrete(observations, model3, algorithmic_parameters)
}

save(r2, r3, file = "~/Dropbox/HyvarinenScore/quick_smcsquare.RData")
load("~/Dropbox/HyvarinenScore/quick_smcsquare.RData")

### diagnostics
get_diagnostics <- function(results){
  diagnostics.df <- data.frame()
  for (r in 1:R){
    accepts <- rep(NA, ncol(observations))
    accepts[results[[r]]$rejuvenation_times] <- results[[r]]$rejuvenation_accepts
    diagnostics.df <- rbind(diagnostics.df, data.frame(time = 1:ncol(observations), ess = results[[r]]$ESS, accepts = accepts, rep = r))
  }
  return(diagnostics.df)
}


diagnostics.df <- get_diagnostics(r2)
g <- ggplot(diagnostics.df, aes(x = time, y = ess, group = rep)) + geom_line() + ylim(0, algorithmic_parameters$Ntheta)
g + geom_vline(xintercept = 3)
matplot(t(observations), type = "l")
diagnostics.df$accepts[is.na(diagnostics.df$accepts)] <- 0
ggplot(diagnostics.df, aes(x = time, y = accepts, group = rep)) + geom_point() + ylim(0.001, 1)
#
# ###  parameters
# thetas.df <- data.frame()
# for (r in 1:5){
#   for (index in 1:ncol(observations)){
#     theta <- data.frame(theta = r2[[r]]$thetas_history[[index]], weight = r2[[r]]$weights_history[[index]], time = index, r = r)
#     names(theta) <- c("sigma", "tau", "r", "weight", "time", "rep")
#     thetas.df <- rbind(thetas.df, theta)
#   }
# }
# thetas.df %>% tail
# ggplot(thetas.df %>% filter(time == 41), aes(x = sigma, weight = weight, group = rep)) + geom_density()
# ggplot(thetas.df %>% filter(time == 41), aes(x = tau, weight = weight, group = rep)) + geom_density()
# ggplot(thetas.df %>% filter(time == 41), aes(x = r, weight = weight, group = rep)) + geom_density()

names(r3[[1]])
r3[[1]]$logevidence

##
modelchoice.df <- data.frame()
for (r in 1:R){
  modelchoice.df <- rbind(modelchoice.df, data.frame(time = 1:ncol(observations),
                                                     logevidence2 = r2[[r]]$logevidence,
                                                     hscore2 = r2[[r]]$hscore,
                                                     logevidence3 = r3[[r]]$logevidence,
                                                     hscore3 = r3[[r]]$hscore, rep = r))
}

modelchoice.df %>% filter(time == 41) %>% summarise(mean2 = mean(logevidence2), sd2 = sd(logevidence2),
                                                    mean3 = mean(logevidence3), sd3 = sd(logevidence3))
modelchoice.df %>% filter(time == 41) %>% summarise(mean2 = mean(hscore2), sd2 = sd(hscore2),
                                                    mean3 = mean(hscore3), sd3 = sd(hscore3))

modelchoice.df %>% group_by(rep) %>% mutate(d = c(0, diff(logevidence2 - logevidence3))) %>% ungroup() %>%
  group_by(time) %>% summarise(s = sd(d))
modelchoice.df %>% group_by(time) %>% summarise(s = sd(logevidence2 - logevidence3))

g <- ggplot(modelchoice.df, aes(x = time, y = logevidence3 - logevidence2, group = rep)) + geom_line(colour = "red")
g + geom_vline(xintercept = 4)

g <- ggplot(modelchoice.df, aes(x = time, y = logevidence2, group = rep)) + geom_line(colour = "orange")
g <- g + geom_line(aes(y = logevidence3), colour = "blue")
g
#
g <- ggplot(modelchoice.df, aes(x = time, y = hscore2 - hscore3, group = rep)) + geom_line(colour = "red")
g + geom_hline(yintercept = 0, linetype = 3)

g <- ggplot(modelchoice.df, aes(x = time, y = hscore2, group = rep)) + geom_line(colour = "orange")
g <- g + geom_line(aes(y = hscore3), colour = "blue")
g


# ggplot(modelchoice.df, aes(x = time, y = hscore, group = rep)) + geom_line()
#sanity check ES
# plot(r3$ESS,type='l', ylim = c(0, algorithmic_parameters$Ntheta))


# ggplot(theta.df, aes(x = sigma, weight = weight)) + geom_histogram()
# ggplot(theta.df, aes(x = tau, weight = weight)) + geom_histogram()
#
# # #sanity check hscore
# # m = min(r1$hscore,r2$hscore,r3$hscore)
# # M = max(r1$hscore,r2$hscore,r3$hscore)
# plot(r3$hscore)
# plot(r3$logevidence)
# # plot(r3$hscore,type='l',ylim = c(m,M))
# # points(r3$hscore,pch=17,col=1,lwd=1)
# # lines(r2$hscore)
# # points(r2$hscore,pch=19,col=1,lwd=1)
# # lines(r1$hscore)
# # points(r1$hscore,pch=15,col=1,lwd=1)
# # legend("bottomleft",legend=c("M3 - Random Walk","M2 - Exponential","M1 - Logistic"),pch=c(17,19,15))
# #
# # #sanity check logevidence
# # m = min(r1$logevidence,r2$logevidence,r3$logevidence)
# # M = max(r1$logevidence,r2$logevidence,r3$logevidence)
# # plot(r3$logevidence,type='l',ylim = c(m,M))
# # points(r3$logevidence,pch=17,col=1,lwd=1)
# # lines(r2$logevidence)
# # points(r2$logevidence,pch=19,col=1,lwd=1)
# # lines(r1$logevidence)
# # points(r1$logevidence,pch=15,col=1,lwd=1)
# # legend("bottomleft",legend=c("M3 - Random Walk","M2 - Exponential","M1 - Logistic"),pch=c(17,19,15))
#
#
#
