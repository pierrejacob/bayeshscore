#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector rtransition_logistic_c(NumericVector Xt, double dt, double sigma, double r, double b){
  RNGScope scope;
  int N = Xt.size();
  NumericVector logXtold = log(Xt);
  int M = dt / 0.001;
  double delta = dt / (double) M;
  double sqrtdelta = sqrt(delta);
  for (int im = 0; im < M; im ++){
    logXtold = logXtold + (r-b*exp(logXtold))*delta + sigma*sqrtdelta*rnorm(N);
  }
  return logXtold;
}
