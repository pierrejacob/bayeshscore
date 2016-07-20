#'@rdname systematic_resampling_n
#'@title Systematic resampling
#'@description Runs a systematic resampling algorithm, taking a uniform variable, normalized weights and a number of desired draws
#'and returning a vector of ancestors
#'@return A vector of ancestors
#'@export

systematic_resampling_n <- function(normalized_weights, N, u){
  return(systematic_resampling_n_(normalized_weights, N, u) + 1)
}