rm(list = ls())
library(HyvarinenSSM)
set.seed(17)
w <- rgamma(16, 4, 2)
w <- w / sum(w)
N <- 24



registerDoMC(cores = 6)
nrep <- 10000
freq_observed <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  a <- systematic_resampling_n(w, N, runif(1))
  freq_observed_ <- tabulate(a, nbins = 16)/24
  freq_observed_
}

summary(abs(apply(freq_observed, 2, mean) - w)/w)

freq_observed <- foreach(irep = 1:nrep, .combine = rbind) %dorng% {
  a <- multinomial_resampling_n(w, N)
  freq_observed_ <- tabulate(a, nbins = 16)/24
  freq_observed_
}

summary(abs(apply(freq_observed, 2, mean) - w)/w)

print("done!")

