#'@rdname multinomial_resampling_n
#'@title Multinomial resampling
#'@description Runs a multinomial resampling algorithm, taking a normalized weights and a number of desired draws
#'and returning a vector of ancestors
#'@return A vector of ancestors
#'@export
multinomial_resampling_n <- function(normalized_weights, N){
  return(multinomial_resampling_n_(normalized_weights, N) + 1)
}
