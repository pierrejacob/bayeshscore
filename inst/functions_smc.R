#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# This file contains the standard particle filter
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# WARNING: the current code only works for univariate states !!!
source("functions_resampling.R")
source("functions_proposal.R")


# This function takes a matrix X (Nx by T) and an ancestral lineage (1 by T)
# and extract the corresponding path
extract_path = function(X,b) {
  T = length(b)
  return (sapply(1:T, FUN = function(i) X[b[i],i]))
}
# This function takes a matrix X (Nx by T) and a matrix of ancestral lineages (Nx by T)
# and extract all the corresponding paths
extract_all_path = function(X,B) {
  N = dim(B)[1]
  t(sapply(1:N,function(i) extract_path(X,B[i,])))
}
# This function runs a particle filter
# it takes a resampling and a proposal functions as additional arguments if needed
# (defautlt is multinomial resampling and kernel transition as proposal)
# The prototype is the following:
# resampling = function(Nx,W[,t-1]) >>> return Nx indexes in 1,...,Nx
# proposal = function(vector_X_old,Y,t,theta) >>> return vector of new particles (Nx by 1)
SMC = function(Nx,y,model,theta,resampling = multinomial,proposal = transition,smoothing = FALSE) {
  T = length(y)
  log_p_y_hat = 0 #initialize estimate of p_theta(y_{1:T}))
  X = matrix(NA,nrow = Nx,ncol = T) #matrix of Nx particles row-wise
  W = matrix(NA,nrow = Nx,ncol = T) #matrix of normalized weights
  A = matrix(NA,nrow = Nx,ncol = T-1) #matrix of ancestors' indexes
  X[,1] = model$rinitial(theta,Nx) #initial step 1
  wt = model$dobs(y[1],X[,1],1,theta)
  log_p_y_hat = log_p_y_hat + mean(wt) #udpate likelihood estimate
  W[,1] = wt/sum(wt) #normalize weights
  #iterate for n = 2, ... T
  for (t in 2:T) {
    A[,t-1] = resampling(Nx,W[,t-1]) #sample the ancestors' indexes
    X[,t] = proposal(X[A[,t-1],t-1],y,t,model,theta) #propagate particles
    wt = model$dobs(y[t],X[,t],t,theta)
    log_p_y_hat = log_p_y_hat + mean(wt) #udpate likelihood estimate  
    W[,t] = wt/sum(wt) #normalize weights
  }
  if (smoothing) {
    B = matrix(NA,nrow = Nx,ncol = T) #initialize matrix of lineages
    B[,T] = 1:Nx
    for (i in 1:Nx) {
      for (t in (T-1):1) {
        B[i,t] = A[B[i,t+1],t]
      }
    }
    return (list(Xsmooth = extract_all_path(X,B),Xfilter = X, W = W,py_hat = exp(log_p_y_hat)))
  }
  else {
    return (list(Xfilter = X, W = W,py_hat = exp(log_p_y_hat)))
  }
}
# Same particle filter with multivariate observations Y
SMCmvtY = function(Nx,y,model,theta,resampling = multinomial,proposal = transition,smoothing = FALSE) {
  T = length(y[1,])
  log_p_y_hat = 0 #initialize estimate of p_theta(y_{1:T}))
  X = matrix(NA,nrow = Nx,ncol = T) #matrix of Nx particles row-wise
  W = matrix(NA,nrow = Nx,ncol = T) #matrix of normalized weights
  A = matrix(NA,nrow = Nx,ncol = T-1) #matrix of ancestors' indexes
  X[,1] = model$rinitial(theta,Nx) #initial step 1
  wt = model$dobs(y[,1],X[,1],1,theta)
  log_p_y_hat = log_p_y_hat + mean(wt) #udpate likelihood estimate
  sum_weight = sum(wt)
  if (sum_weight>0) {
    W[,1] = wt/sum_weight #normalize weights
  }
  else {
    W[,1] = rep(1/Nx,Nx) #WARNING:  NOT SURE ABOUT THIS CASE !!! IS IT DUE TO UNDERFLOW ??
  }
  #iterate for n = 2, ... T
  for (t in 2:T) {
    A[,t-1] = resampling(Nx,W[,t-1]) #sample the ancestors' indexes
    for (n in 1:Nx) {
      X[n,t] = proposal(X[A[n,t-1],t-1],y,t,model,theta) #propagate particles
    }
    wt = model$dobs(y[,t],X[,t],t,theta)
    log_p_y_hat = log_p_y_hat + mean(wt) #udpate likelihood estimate  
    sum_weight = sum(wt)
    if (sum_weight>0) {
      W[,t] = wt/sum_weight #normalize weights
    }
    else {
      W[,t] = rep(1/Nx,Nx) #WARNING:  NOT SURE ABOUT THIS CASE !!! IS IT DUE TO UNDERFLOW ??
    }
  }
  if (smoothing) {
    B = matrix(NA,nrow = Nx,ncol = T) #initialize matrix of lineages
    B[,T] = 1:Nx
    for (i in 1:Nx) {
      for (t in (T-1):1) {
        B[i,t] = A[B[i,t+1],t]
      }
    }
    return (list(Xsmooth = extract_all_path(X,B),Xfilter = X, W = W,py_hat = exp(log_p_y_hat)))
  }
  else {
    return (list(Xfilter = X, W = W,py_hat = exp(log_p_y_hat)))
  }
}

