rm(list = ls())
library(HyvarinenSSM)
set.seed(17)

rtransition = function(Xt,t,theta){
  sigma = theta[1]
  r = theta[3]
  b = theta[4]
  N = nrow(Xt)
  logXtold = log(Xt)
  dt = data_kangaroo["time",t] - data_kangaroo["time",t-1]
  M = dt/0.001
  delta = dt/M #Euler method discretization step size
  sqrtdelta <- sqrt(delta)
  for (i in 1:M) {
    logXtold = logXtold + (r-b*exp(logXtold))*delta + sigma*sqrtdelta*rnorm(N)
    # logXtold = logXtnew
  }
  return (matrix(exp(logXtold),nrow = N))
}

library(Rcpp)

cppFunction("
NumericVector rtransition_rcpp(NumericVector Xt, double dt, double sigma, double r, double b){
  int N = Xt.size();
  NumericVector logXtold = log(Xt);
  int M = dt / 0.001;
  double delta = dt / (double) M;
  double sqrtdelta = sqrt(delta);
  for (int im = 0; im < M; im ++){
    logXtold = logXtold + (r-b*exp(logXtold))*delta + sigma*sqrtdelta*rnorm(N);
  }
  return logXtold;
}
")


nobservations <- 20
model <- get_model_kangarooLogistic()
#simulate data
sim = simulateData(model, theta = c(0.1,0.1,0.1,1), nobservations)
X = sim$X
Y = sim$Y
theta = c(0.1,0.1,0.1,1)
Xt <- model$rinitial(theta = theta, N = 1024)
rtransition(Xt = Xt, t = 2, theta = theta)
t <- 2
rtransition_rcpp(Xt[,1], data_kangaroo["time",t] - data_kangaroo["time",t-1], theta[1], theta[3], theta[4])
log(Xt[,1])


library(microbenchmark)
microbenchmark(  rtransition(Xt = Xt, t = 2, theta = theta),
                 rtransition_logistic(Xt, 2, theta),
                 times = 200)
