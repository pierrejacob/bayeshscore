#------------------------------------------------------------------------------#
#----------------------  Some resampling algorithms ---------------------------#
#------------------------------------------------------------------------------#
#'@rdname multinomial_resampling_n
#'@title Multinomial resampling
#'@description Runs a multinomial resampling algorithm, taking a normalized weights and a number of desired draws
#'and returning a vector of ancestors
#'@return A vector of ancestors
#'@export
multinomial_resampling_n <- function(normalized_weights, N){
  return(multinomial_resampling_n_(normalized_weights, N) + 1)
}
#------------------------------------------------------------------------------#
#'@rdname systematic_resampling_n
#'@title Systematic resampling
#'@description Runs a systematic resampling algorithm, taking a vector of normalized weights, a number of desired draws, and a uniform variable
#'and returning a vector of ancestors
#'@return A vector of ancestors
#'@export
systematic_resampling_n <- function(normalized_weights, N, u){
  return(systematic_resampling_n_(normalized_weights, N, u) + 1)
}
#------------------------------------------------------------------------------#
#'@rdname ssp_resampling_n
#'@title SSP resampling
#'@description Runs an SSP resampling algorithm, taking a vector of normalized weights (length N) and a vector of uniform variables (length N)
#'and returning a vector of ancestors (length N)
#'@return A vector of ancestors
#'@export
ssp_resampling_n <- function(normalized_weights, u_vector){
  return(SSP_resampling_n_(normalized_weights, u_vector) + 1)
}
