#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# This file contains some generic functions for models
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

# This function simulates data (T titmesteps) from a given SSM model with parameter theta
simulateData = function(model,theta,T) {
  if (is.null(model$robs)) {
    print("no available simulator for observations")
    return (NULL)
  }
  else {
    X = vector("numeric",T)
    Y = vector("numeric",T)
    X[1] = model$rinitial(theta,1)
    Y[1] = model$robs(X[1],1,theta)
    for (t in 2:T) {
      X[t] = model$rtransition(X[t-1],t,theta)
      Y[t] = model$robs(X[t],t,theta)
    }
    return (list(X = X, Y = Y))
  }
}

simulateDataKangaroo = function(kangaroomodel,theta,timesteps) {
  if (is.null(kangaroomodel$robs)) {
    print("no available simulator for observations")
    return (NULL)
  }
  else {
    T = length(timestep)
    X = vector("numeric",T)
    Y = array(NA,dim = c(2,T))
    X[1] = kangaroomodel$rinitial(theta,1)
    Y[,1] = kangaroomodel$robs(X[1],1,theta)
    for (t in 2:T) {
      dt = timestep[t]-timestep[t-1]
      X[t] = kangaroomodel$rtransition(X[t-1],dt,theta)
      Y[,t] = kangaroomodel$robs(X[t],t,theta)
    }
    return (list(X = X, Y = Y))
  }
}