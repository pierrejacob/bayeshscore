#-----------------------------------------------------------------------------------#
# Functions to estimate the derivatives of the density by following
# the approach of Sasaki, Noh,  Sugiyama (2015).
# We limit ourselved to the univariate case for simplicity.
# We use gaussian RBF to model p'(y) and p''(y).
#-----------------------------------------------------------------------------------#
#'@rdname get_estimate_RBF
#'@title get_estimate_RBF
#'@description Evaluate the RBF model at y given the fitted parameters and the training set ys.
#'@export
get_estimate_RBF = function(y, ys, theta, sigma2){
  return (sapply(y,function(x)sum(theta*exp(-((x - ys)^2)/(2*sigma2)))))
}

#'@rdname get_CV_RBF
#'@title get_CV_RBF
#'@description Finds the optimal lambda and sigma2 for p'(y) or p''(y) via cross-validation.
#'The function assumes length(ys) is a multiple of kfold.
#'@export
get_CV_RBF = function(ys, kfold, range_lambda, range_sigma2, order) {
  Ntraining = length(ys)/kfold
  cv_loss = Inf
  cv_lambda = range_lambda[1]
  cv_sigma2 = range_sigma2[1]
  if (order == 1){loss_function = loss_d1; get_h = get_h1_RBF_cpp}
  if (order == 2){loss_function = loss_d2; get_h = get_h2_RBF_cpp}
  for (lambda in range_lambda){
    for (sigma2 in range_sigma2){
      loss = 0
      for (k in 1:kfold) {
        training_indexes = ((k-1)*Ntraining+1):(k*Ntraining)
        ytrain = ys[training_indexes]
        ytest = ys[-training_indexes]
        G = get_G_RBF_cpp(ytrain, sigma2)
        h = get_h(ytrain, sigma2)
        theta = solve(G + diag(lambda, Ntraining, Ntraining), -h)
        loss = loss + loss_function(theta, sigma2, G, ytrain, ytest)
      }
      loss = loss/kfold
      if (loss < cv_loss){
        cv_loss = loss
        cv_lambda = lambda
        cv_sigma2 = sigma2
      }
    }
  }
  return (list(sigma2 = cv_sigma2, lambda = cv_lambda))
}

#'@rdname get_derivative_RBF
#'@title get_derivative_RBF
#'@description Estimates p'(y) or p''(y) via RBF regression.
#'The function assumes length(ys) is a multiple of kfold.
#'@export
get_derivative_RBF = function(ys, lambda, sigma2, order, ynew = NULL) {
  if (order == 1){get_h = get_h1_RBF_cpp}
  if (order == 2){get_h = get_h2_RBF_cpp}
  N = length(ys)
  G = get_G_RBF_cpp(ys, sigma2)
  h = get_h(ys, sigma2)
  theta = solve(G + diag(lambda, N, N), -h)
  estimate = NULL
  if (!is.null(ynew)){
    estimate = get_estimate_RBF(ynew, ys, theta, sigma2)
  }
  return (list(theta = theta, sigma2 = sigma2, lambda = lambda, estimate = estimate))
}

#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
# Variation with locally weighted loss
#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
#'@rdname get_derivative_RBFlocal
#'@title get_derivative_RBFlocal
#'@description Estimates derivative of p(y) up to order 2 by fitting a constant via a local loss.
#'@export
get_derivative_RBFlocal = function(ys, ystar, sigma2star, order) {
  return (get_derivative_cpp(ys, ystar, sigma2star, order))
}
