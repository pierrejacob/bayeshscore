#'@export
# This function computes the ESS given a set of normalized weights
getESS = function(normalized_weights) {
  return (1/sum(normalized_weights^2))
}
