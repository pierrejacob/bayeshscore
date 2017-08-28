


n = 100
A = matrix(rnorm(n*n),ncol=n)
b = matrix(rnorm(n),ncol=1)

lambdas = seq(0,1,0.01)

f1 = function(){
  for (i in 1:length(lambdas)){
    solve(A+diag(lambdas[i],n),b)
  }
}

f2 = function(){
  svdA = svd(A)
  for (i in 1:length(lambdas)){
    svdA$v %*% diag(1/(svdA$d+lambdas[i])) %*% t(svdA$u) %*% b
  }
}

microbenchmark::microbenchmark(f1(),f2())
