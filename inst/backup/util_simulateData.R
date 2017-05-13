#'@rdname simulateData
#'@title simulateData
#'@description This function generates artificial data from a model. 'additional_parameters' is a list of n lists, where n is the number of observations that are generated (examples of additional parameters include step size for simulating SDEs)
#'@export
simulateData = function(model,theta,nobservations) {
  X <- matrix(nrow = model$dimX, ncol = nobservations)
  Y <- matrix(nrow = model$dimY, ncol = nobservations)
  X[,1] <- model$rinitial(theta,1)
  Y[,1] <- model$robs(matrix(X[,1], nrow = 1),1,theta)
  for (t in 2:nobservations) {
    X[,t] = model$rtransition(X[,t-1,drop=FALSE], t, theta)
    Y[,t] = model$robs(X[,t-1,drop=FALSE], t, theta)
  }
  return (list(X = X, Y = Y))
}

