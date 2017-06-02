#'@export
#'@rdname get_derivatives_from_genD
#'@title get_derivatives_from_genD
#'@description This function extracts the jacobian and the diagonal of the hessians from the output of genD (cf. numDeriv)
#'@export
get_derivatives_from_genD <- function(genD_output,dimY){
  indexesofdiagonalcoeff = rep(1+dimY,dimY)
  if (!is.null(dimY)) {
    if (dimY > 1){
      increment = 2
      for (i in 2:dimY){
        indexesofdiagonalcoeff[i] = indexesofdiagonalcoeff[i-1] + increment
        increment = increment + 1
      }
    }
  }
  return(list(jacobian = genD_output[,1:dimY,drop=FALSE], hessiandiag = genD_output[,indexesofdiagonalcoeff,drop=FALSE]))
}
