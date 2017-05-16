#------------------------------------------------------------------------------#
#----------------------  Some useful functions for tempering ------------------#
#------------------------------------------------------------------------------#
#'@export
#'@rdname search_gamma
#'@title search_gamma
#'@description This function finds the optimal gamma for tempering given some objective function f (e.g. ESS)
#' WARNING: f has to be a continuous function
#'@export
search_gamma <- function(current_gamma, f, objective, maxsteps = 1000, tolerance = 1e-2){
  if ((f(current_gamma) < objective)|| f(1) > objective){
    print("problem! there's no solution to the binary search")
  }
  attempt <- 1
  current_size <- (1 - current_gamma)/2
  fattempt <- f(attempt)
  istep <- 0
  while (!(fattempt >= objective && fattempt < objective+tolerance) && (istep < maxsteps)){
    istep <- istep + 1
    if (fattempt > objective){
      attempt <- attempt + current_size
      fattempt <- f(attempt)
      current_size <- current_size / 2
    } else {
      attempt <- attempt - current_size
      fattempt <- f(attempt)
      current_size <- current_size / 2
    }
  }
  return(list(x = attempt, f = fattempt))
}
#------------------------------------------------------------------------------#
#'@export
#'@rdname getESS
#'@title getESS
#'@description This function computes the ESS given a set of normalized weights
#'@export
getESS = function(normalized_weights) {
  return (1/sum(normalized_weights^2))
}
