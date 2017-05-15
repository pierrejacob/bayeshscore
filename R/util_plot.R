#------------------------------------------------------------------------------#
#---------------------------  Some useful for plotting graphs -----------------#
#------------------------------------------------------------------------------#
#'@rdname simulateData
#'@title simulateData
#'@description This function generates artificial data from a model.
#'@export
simulateData = function(model,theta,nobservations) {
  X <- matrix(nrow = model$dimX, ncol = nobservations)
  Y <- matrix(nrow = model$dimY, ncol = nobservations)
  X[,1] <- model$rinitial(theta,1)
  Y[,1] <- model$robs(matrix(X[,1], nrow = 1),1,theta)
  for (t in 2:nobservations) {
    X[,t] = model$rtransition(X[,t-1,drop=FALSE], t, theta)
    Y[,t] = model$robs(X[,t,drop=FALSE], t, theta)
  }
  return (list(X = X, Y = Y))
}
#------------------------------------------------------------------------------#
#'@rdname setmytheme
#'@title Customize graphical settings
#'@description This function customizes the theme used by ggplot.
#'@export
setmytheme <- function(){
  theme_set(theme_bw())
  theme_update(axis.text.x = element_text(size = 20),
               axis.text.y = element_text(size = 20),
               axis.title.x = element_text(size = 25, margin=margin(20,0,0,0)),
               axis.title.y = element_text(size = 25, angle = 90, margin = margin(0,20,0,0)),
               legend.text = element_text(size = 20),
               legend.title = element_text(size = 20),
               title = element_text(size = 30),
               strip.text = element_text(size = 25),
               strip.background = element_rect(fill="white"),
               panel.margin = unit(2, "lines"),
               legend.position = "bottom")
}
