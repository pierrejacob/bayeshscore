#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector rtransition_logistic_c(NumericVector Xt, double delta_t, double dt, double sigma, double r, double b){
  RNGScope scope;
  int N = Xt.size();
  NumericVector logXtold = log(Xt);
  int M = delta_t / dt;
  double sqrtdt = sqrt(dt);
  for (int im = 0; im < M; im ++){
    logXtold = logXtold + (r-b*exp(logXtold))*dt + sigma*sqrtdt*rnorm(N);
  }
  return exp(logXtold);
}
