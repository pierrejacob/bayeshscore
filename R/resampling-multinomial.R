#' @export

multinomial_resampling_n <- function(normalized_weights, N){
  return(multinomial_resampling_n_(normalized_weights, N) + 1)
}