#'@export
#'@rdname repeat_column
#'@title repeat_column
#'@description This function repeats n times a column vector X to form a nrow(X) by n matrix
#'@export
repeat_column<-function(n,x){
  matrix(rep(x,n), ncol=n)
}
