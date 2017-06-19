library(Rcpp)

N = 1000
test = rnorm(N)

fr = function(x){
  return (x %o% x)
}

cppFunction('NumericMatrix fc(NumericVector x) {
  int N = x.length();
  NumericMatrix result(N);
  for (int i = 0; i < N; i ++){
    for (int j = 0; j <= i; j ++){
      result(i,j) = x(i)*x(j);
    }
  }
  for (int i = 0; i < N; i ++){
    for (int j = (i+1); j < N; j ++){
      result(i,j) = result(j,i);
    }
  }
  return result;
}')


microbenchmark::microbenchmark(fr(test),fc(test))
