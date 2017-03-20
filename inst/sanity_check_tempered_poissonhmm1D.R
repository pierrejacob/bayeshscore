rm(list = ls())
library(doParallel)
library(HyvarinenSSM)
library(gridExtra)
library(numDeriv)
set.seed(19)
module_tree <<- Module("module_tree", PACKAGE = "HyvarinenSSM")
TreeClass <<- module_tree$Tree

# Define model and data
nobservations <- 5
model <- get_model_poissonhmm()
true_theta = 0.29
sim = simulateData(model, true_theta, nobservations)
X = sim$X
Y = sim$Y
observations <- matrix(Y, nrow = model$dimY)# observations in a matrix of dimensions dimy x nobservations

# Plot data
observations.df = data.frame(time = 1:nobservations, X = t(X),Y = t(Y))
g = ggplot(observations.df, aes(x = time)) +
  geom_point(aes(y=Y),size=2) +
  geom_line(aes(y=X),linetype=2) +
  xlab("\n Time")
plot(g)

# Define algorithmic parameters for each model
Ntheta = 2^10
algorithmic_parameters = list(Ntheta = Ntheta,
                              resampling = function(normw) systematic_resampling_n(normw, length(normw), runif(1)),
                              progress = TRUE, nmoves = 1)
algorithmic_parameters$Nx = NULL # this uses adaptive Nx, starting with Nx = 128
algorithmic_parameters$min_acceptance_rate = 0.45
algorithmic_parameters$ess_threshold = 0.5

repl = 5
results.df = data.frame()
posterior.df = data.frame()
results_temp.df = data.frame()
posterior_temp.df = data.frame()
for (r in 1:repl){
  # Tempered SMC^2
  results_temp = hscore_discrete(observations, model, algorithmic_parameters)
  # Regular SMC^2
  results = hscore_discrete_no_tempering(observations, model, algorithmic_parameters)
  results.df = rbind(results.df, data.frame(time = 1:ncol(observations),
                                            logevidence = results$logevidence,
                                            hscore = results$hscore,
                                            ess = results$ESS,
                                            rep = r,
                                            sim = 1))
  results_temp.df = rbind(results_temp.df, data.frame(time = 1:ncol(observations),
                                                      logevidence = results_temp$logevidence,
                                                      hscore = results_temp$Hscore,
                                                      rep = r,
                                                      sim = 1))
  posterior.df = rbind(posterior.df, data.frame(theta = results$thetas[,1],
                                                w = results$thetanormw,
                                                rep = r,
                                                sim = 1))
  posterior_temp.df = rbind(posterior_temp.df, data.frame(theta = results_temp$thetas_history[[nobservations+1]][,1],
                                                          w = results_temp$normw_history[[nobservations+1]],
                                                          rep = r,
                                                          sim = 1))
}


#Compute exact likelihood
l = model$lambda
lx = model$lambdaX
p_Y_1_t_theta = function(Y_1_t,theta){
  t = ncol(Y_1_t)
  all_paths = expand.grid(rep(list(0:1), t))
  if (t==1){
    temp = 0
    for (i in 1:(2^t)){
      X = all_paths[i,1:t]
      rate = l+X*lx
      temp = temp + prod(exp(-(rate))*(rate^Y_1_t)/factorial(Y_1_t))*(1/2) + 0*theta
      #the "0*theta" term (artificially) allows integrate to recognize this as a function of theta
    }
    return (temp)
  }
  else {
    temp = 0
    for (i in 1:(2^t)){
      X = all_paths[i,1:t]
      rate = l+X*lx
      temp = temp + prod(exp(-(rate))*(rate^Y_1_t)/factorial(Y_1_t))*(1/2)*(theta^sum((1-X[2:t])*(1-X[1:(t-1)])+X[2:t]*X[1:(t-1)]))*((1-theta)^sum(X[2:t]*(1-X[1:(t-1)])+(1-X[2:t])*X[1:(t-1)]))
    }
    return (temp)
  }
}
#Compute exact evidence
p_Y_1_t = function(Y_1_t) {
  likelihood = function(theta) p_Y_1_t_theta(Y_1_t,theta)
  return (integrate(likelihood,0,1)$value)
}
#Compute exact log-evidence
logevidence_exact = rep(NA,nobservations)
for (t in 1:nobservations) {
  logevidence_exact[t] = log(p_Y_1_t(Y[,1:t,drop=FALSE]))
}
#Compute exact posterior
normalizing = integrate(function(theta) p_Y_1_t_theta(Y,theta),0,1)$value
post = function(theta){
  return ((p_Y_1_t_theta(Y,theta)/normalizing)*(theta<=1)*(theta>=0))
}
#Compute exact predictive
py_t_t_1_func = list()
for (t in 1:nobservations){
  py_t_t_1_func[[t]] = function(Yt){
    if (t==1){
      return (p_Y_1_t(matrix(Yt,ncol=1)))
    }
    else{
      Y_1_t = Y[,1:t,drop=FALSE]
      Y_1_t[,t] = Yt
      return (p_Y_1_t(Y_1_t)/p_Y_1_t(Y_1_t[,1:(t-1),drop=FALSE]))
    }
  }
}
hincrement = function(k,a,b,d,y,t) {
  ek = rep(0,d)
  ek[k] = 1
  q = py_t_t_1_func[[t]]
  if (y[k]==b[k]) {
    qy = q(y)
    qy_minusek = q(y-ek)
    return (-2*(qy-qy_minusek)/qy_minusek)
  }
  else {
    if (y[k]==a[k]) {
      qy = q(y)
      qy_plusek = q(y+ek)
      return (2*(qy_plusek-qy)/qy + ((qy_plusek-qy)/qy)^2)
    }
    else {
      qy = q(y)
      qy_minusek = q(y-ek)
      qy_plusek = q(y+ek)
      return (2*((qy_plusek-qy)/qy-(qy-qy_minusek)/qy_minusek) + ((qy_plusek-qy)/qy)^2)
    }
  }
}
#Compute exact hscore
exact_hscore = rep(NA,nobservations)
for (t in 1:nobservations){
  for (k in 1:model$dimY){
    exact_hscore[t] = hincrement(k,model$lower,model$upper,model$dimY,Y[,t],t)
  }
}
exact_preq_hscore = cumsum(exact_hscore)


#plot Posterior
g11 <- ggplot(posterior.df, aes(x = theta, weight = w, group = rep)) +
  stat_function(fun = function(theta)post(theta), colour = "red",size=1.5,linetype=1) +
  geom_density(color = "blue",linetype="dashed") + ylab("") +
  geom_density(data = posterior_temp.df, aes(x = theta, weight = w, group = rep),color = "blue")
grid.arrange(g11,ncol = 1, nrow = 1)


#Append exact results to dataframe
results.df = rbind(results.df, data.frame(time = 1:nobservations,
                                          logevidence = logevidence_exact,
                                          hscore = exact_preq_hscore,
                                          ess = rep(1,nobservations),
                                          rep = 0,
                                          sim = 0))

#CHECK LOG-EVIDENCE (red is exact)
g <- ggplot() +
  geom_point(aes(x = time, y = logevidence,color="Exact"),data = subset(results.df,sim==0), colour = "red",size=4) +
  geom_line(aes(x = time, y = logevidence),data = subset(results.df,sim==1), colour="blue",size = 2,linetype="dashed") +
  geom_line(aes(x = time, y = logevidence, group = rep),data = results_temp.df, size = 1, color = "blue")
plot(g)

#CHECK HSCORE (red is exact)
g <- ggplot() +
  geom_point(aes(x = time, y = hscore,color="Exact"),data = subset(results.df,sim==0), colour = "red",size=4) +
  geom_line(aes(x = time, y = hscore),data = subset(results.df,sim==1), colour="blue",size = 2,linetype="dashed") +
  geom_line(aes(x = time, y = hscore, group = rep),data = results_temp.df, size = 1, color = "blue")
plot(g)
