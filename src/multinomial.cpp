#include <RcppEigen.h>
#include "multinomial.h"
using namespace Rcpp;
using namespace std;

inline int randWrapper(const int n) {
  RNGScope scope;
  NumericVector ru = runif(1);
  return floor(ru(0)*n);
}

// [[Rcpp::export]]
IntegerVector multinomial_resampling_n_(const NumericVector & weights, int ndraws){
  RNGScope scope;
  int nparticles = weights.size();
  IntegerVector ancestors(ndraws);
  NumericVector cumsumw = cumsum(weights);
  NumericVector uniforms = runif(ndraws);
  double sumw = cumsumw(nparticles - 1);
  double lnMax = 0;
  int j = nparticles;
  for (int i = ndraws; i > 0; i--){
    lnMax += log(uniforms(i-1)) / i;
    uniforms(i-1) = sumw * exp(lnMax);
    while (j > 0 && uniforms(i-1) < cumsumw(j-1)){
      j --;
    }
    ancestors(i-1) = j;
  }
  std::random_shuffle(ancestors.begin(), ancestors.end(), randWrapper);
  return ancestors;
}

