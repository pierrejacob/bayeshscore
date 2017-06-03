#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;

// Simulate latent state transition of model 1 for kangaroo counts (cf. Knape and Valpine, 2012)
// [[Rcpp::export]]
NumericMatrix rtransition_logistic_cpp(const NumericMatrix& Xold, const double& delta_t, const double& dt,
                                       const double& sigma, const double& r, const double& b){
  RNGScope scope;
  int N = Xold.ncol();
  NumericMatrix Xnew(1,N);
  Xnew(0,_) = log(Xold(0,_));
  int M = delta_t / dt;
  double sqrtdt = sqrt(dt);
  for (int im = 0; im < M; im ++){
    Xnew(0,_) = Xnew(0,_) + (r-b*exp(Xnew(0,_)))*dt + sigma*sqrtdt*rnorm(N);
  }
  Xnew(0,_) = exp(Xnew(0,_));
  return Xnew;
}

