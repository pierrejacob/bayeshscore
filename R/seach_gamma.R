# f has to be a decreasing continuous function
#'@export
seach_gamma <- function(current_gamma, f, objective, maxsteps = 1000, tolerance = 1e-2){
  if ((f(current_gamma) < objective)|| f(1) > objective){
    print("problem! there's no solution to the binary search")
  }
  attempt <- 1
  current_size <- (1 - current_gamma)/2
  fattempt <- f(attempt)
  istep <- 0
  while (!(fattempt > objective-tolerance && fattempt < objective+tolerance) && (istep < maxsteps)){
    istep <- istep + 1
    if (fattempt > objective-tolerance){
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
