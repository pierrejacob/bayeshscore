library(HyvarinenSSM)

f <- get_independent_normal_proposal()
library(mvtnorm)

n <- 1e6
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

colMeans(thetas)
mean_theta

cov(thetas)
cov_theta


normw <- rep(1/n, n)
prop <- f(t(thetas), normw)

#
covariance = cov.wt(thetas, wt = normw, method = "ML")
mean_t = covariance$center
cov_t = covariance$cov + diag(rep(10^(-4)/ncol(thetas)), ncol(thetas)) # increased a bit the diagonal to prevent degeneracy effects)
#

X <- prop$r(1e6)
rowMeans(X)
mean_theta
cov(t(X))
cov_theta

evals <- prop$d(X, log = T)
evals_mvtnorm <- dmvnorm(t(X), mean = mean_t, sigma = cov_t, log = TRUE)
summary(abs(evals - evals_mvtnorm))
