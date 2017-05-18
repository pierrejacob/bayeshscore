library(HyvarinenSSM)

f <- get_mixture_proposal(nclust = 5)
library(mvtnorm)

n <- 1e3
d <- 5

mean_theta <- rnorm(d)
cov_theta <- diag(1, d, d)
for (i in 1:d){
  for (j in 1:d){
    cov_theta[i,j] <- 0.8^(abs(i-j))
  }
}
cov_theta
thetas <- rmvnorm(n, mean_theta, cov_theta)
thetas <- t(thetas)

normw <- rep(1/n, n)

prop <- f(thetas, normw)

X <- prop$r(Ntheta = 1e6)
hist(X[1,], prob = TRUE, nclass = 100, xlim = c(-6,4), ylim = c(0, 1))
hist(thetas[1,], prob = T, nclass = 30, add = TRUE, xlim = c(-6,4), col = rgb(1,0,0,alpha = 0.4))

hist(X[3,], prob = TRUE, nclass = 100, xlim = c(-6,4), ylim = c(0, 1))
hist(thetas[3,], prob = T, nclass = 30, add = TRUE, xlim = c(-6,4), col = rgb(1,0,0,alpha = 0.4))


rowMeans(X)
mean_theta
cov(t(X))
cov_theta

