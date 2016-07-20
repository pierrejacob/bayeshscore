#'@export
simulateData = function(model,theta,nobservations) {
  X <- matrix(nrow = model$dimX, ncol = nobservations)
  Y <- matrix(nrow = model$dimY, ncol = nobservations)
  X[,1] <- model$rinitial(theta,1)
  Y[,1] <- model$robs(matrix(X[,1], nrow = 1),1,theta)
  for (t in 2:nobservations) {
    X[,t] = model$rtransition(matrix(X[,t-1], nrow = 1), t, theta)
    Y[,t] = model$robs(matrix(X[,t], nrow = 1), t, theta)
  }
  return (list(X = X, Y = Y))
}
