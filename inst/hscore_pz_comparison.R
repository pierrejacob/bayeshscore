rm(list = ls())
library(HyvarinenSSM)
library(doMC)

set.seed(17)

model.pz <- get_model_pz4()
model.pz$theta = c(0.7,0.5,0.25,0.3)
nobservations <- 100
datafile <- paste0("~/Dropbox/Harvard/Grant/NSF DMS Fall 2016/pzdata.T", nobservations, ".RData")
## generate data
# sim = simulateData(model.pz, theta = c(0.7,0.5,0.25,0.3), nobservations)
# save(sim, file = datafile)
## load data
load(datafile)
X = sim$X
Y = sim$Y

observations <- matrix(Y, nrow = model.pz$dimY)# observations in a matrix of dimensions dimy x nobservations

model.pz5 <- get_model_pz5()

algorithmic_parameters <- list(Ntheta = 2^8, Nx = 2^8, nmoves = 3,
                               resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
                               progress = TRUE)

nrep <- 4

load("~/Dropbox/Harvard/Grant/NSF DMS Fall 2016/hscore.pz4.T100Nt1024Nx512.RData")
load("~/Dropbox/Harvard/Grant/NSF DMS Fall 2016/hscore.pz5.T100Nt1024Nx512.RData")
load("~/Dropbox/Harvard/Grant/NSF DMS Fall 2016/hscore.pz6.T100Nt1024Nx512.RData")

scores.df <- foreach (irep = 1:nrep, .combine = rbind) %dopar% {
  hscores <- c(results.pz4[[irep]]$hscore, results.pz5[[irep]]$hscore, results.pz6[[irep]]$hscore)
  logscores <- -c(results.pz4[[irep]]$logevidence, results.pz5[[irep]]$logevidence, results.pz6[[irep]]$logevidence)
  models <- rep(c("PZ4", "PZ5", "PZ6"), each = nobservations)
  time <- rep(1:nobservations, 3)
  i <- rep(irep, 3*nobservations)
  data.frame(time = time, irep = i, model = models, Hscore = hscores, logscore = logscores)
}

scores.df %>% tail

ggplot(scores.df, aes(x = time, y = Hscore, group = interaction(irep, model), colour = model)) + geom_line()
ggplot(scores.df, aes(x = time, y = logscore, group = interaction(irep, model), colour = model)) + geom_line()

scores5versus6.df <- foreach (irep = 1:nrep, .combine = rbind) %dopar% {
  hscores5 <- results.pz5[[irep]]$hscore
  hscores6 <- results.pz6[[irep]]$hscore
  hscores5minus6 <- hscores5 - hscores6
  logscores5 <- -results.pz5[[irep]]$logevidence
  logscores6 <- -results.pz6[[irep]]$logevidence
  logscores5minus6 <- logscores5 - logscores6
  time <- 1:nobservations
  i <- rep(irep, nobservations)
  data.frame(time = time, irep = i, hscores5 = hscores5, hscores6 = hscores6, hscores5minus6 = hscores5minus6,
             logscores5 = logscores5, logscores6 = logscores6, logscores5minus6 = logscores5minus6)
}
scores5versus6.df %>% head
g <- ggplot(scores5versus6.df %>% group_by(time) %>% summarise(Hcomparison = mean(hscores5minus6)),
       aes(x = time, y = Hcomparison)) + geom_line()
g <- g + ylab("H-score difference")
g
ggsave(filename = "~/Dropbox/Harvard/Grant/NSF DMS Fall 2016/pzHscore.pdf", plot = g, height = 5, width = 10)

g <- ggplot(scores5versus6.df %>% group_by(time) %>% summarise(logcomparison = mean(logscores5minus6)),
       aes(x = time, y = logcomparison)) + geom_line()
g <- g + ylab("log-score difference")
g
ggsave(filename = "~/Dropbox/Harvard/Grant/NSF DMS Fall 2016/pzlogscore.pdf", plot = g, height = 5, width = 10)


# ggplot(scores5versus6.df %>% filter(irep <= 1), aes(x = time, y = hscores5minus6, group = irep)) + geom_line()
# ggplot(scores5versus6.df %>% filter(irep <= 1), aes(x = time, y = logscores5minus6, group = irep)) + geom_line()

# Hdf <- melt(hscores4- hscores5)
# names(Hdf) <- c("time", "rep", "Hscore")
# logdf <- melt(evidences5- evidences4)
# names(logdf) <- c("time", "rep", "logscore")
# score.df <- merge(Hdf, logdf, by = c("time", "rep"))
# CoupledCPF::setmytheme()
# g <- ggplot(score.df, aes(x = time, y = Hscore, group = rep)) + geom_line() + ylab("H-score")
# g
# # ggsave(filename = "~/Dropbox/Harvard/Grant/NSF DMS Fall 2016/pzHscore.pdf", plot = g, height = 5, width = 10)
# g <- ggplot(score.df, aes(x = time, y = logscore, group = rep)) + geom_line() + ylab("logarithmic score")
# g
# # ggsave(filename = "~/Dropbox/Harvard/Grant/NSF DMS Fall 2016/pzlogscore.pdf", plot = g, height = 5, width = 10)
#
# # matplot(evidences, type = "l")
# # matplot(ESSs, type = "l")
#
# # length(unique(results[[1]]$thetas))
# #sanity check
# # par(mfrow=c(3,1))
# # plot(result$ESS,type='l')
# # plot(result$logevidence,type='l')
# # plot(result$hscore,type='l')
# # par(mfrow=c(1,1))
#
# # print(result$logevidence)
# # print(result$hscore)
