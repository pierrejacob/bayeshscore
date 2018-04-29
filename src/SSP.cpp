#include <RcppEigen.h>
#include "SSP.h"
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerVector SSP_resampling_n_(const NumericVector & weights, const NumericVector & u, const double tol = 1e-15){
  RNGScope scope;
  int N = weights.size();
  int n,m,k;
  NumericVector Y(N);
  IntegerVector nb_offsprings(N);
  IntegerVector ancestors(N);
  if (N==1){
    ancestors(0) = 0;
  } else {
    for(int i = 0; i < N; i++){
      Y(i) = N*weights(i);
      nb_offsprings(i) = floor(Y(i));
    }
    n = 0;
    m = 1;
    for(k = 0; k < N; k++){
      if ((m+1) > N){
        break;
      }
      double delta_n = nb_offsprings(n)+1-Y(n);
      double delta_m = Y(m)-nb_offsprings(m);
      double epsilon_n = Y(n)-nb_offsprings(n);
      double epsilon_m = nb_offsprings(m)+1-Y(m);
      //// Test if Y(n) or Y(m) is an integer
      bool Yn_is_integer = ((delta_n < tol) || (epsilon_n < tol));
      bool Ym_is_integer = ((delta_m < tol) || (epsilon_m < tol));
      //// Deal with easy cases first: either Y(n) or Y(m) (or both) is already an integer
      if (Yn_is_integer && Ym_is_integer){
        n = m+1;
        m = m+2;
      } else if (Yn_is_integer){
        n = m;
        m = m+1;
      } else if (Ym_is_integer){
        m = m+1;
      } else {
        double delta = min(delta_n, delta_m);
        double epsilon = min(epsilon_n, epsilon_m);
        //// Keeping track of delta_n, delta_m, epsilon_n, epsilon_m prevents rounding errors
        if (u(k) < epsilon/(epsilon+delta)){
          if (std::abs(delta_n-delta_m)<tol){
            // Test delta_n == delta_m while accounting for numerical rounding imprecisions
            nb_offsprings(n) += 1;
            n = m+1;
            m = m+2;
          } else if (delta_n < delta_m){
            Y(m) = Y(m) - delta;
            nb_offsprings(n) += 1;
            n = m;
            m = m+1;
          } else if (delta_n > delta_m){
            Y(n) = Y(n) + delta;
            m = m+1;
            // Note that nb_offsprings(m) is already set at the correct value
          }
        } else {
          if (std::abs(epsilon_n-epsilon_m)<tol){
            // Test epsilon_n == epsilon_m while accounting for numerical rounding imprecisions
            nb_offsprings(m) += 1;
            n = m+1;
            m = m+2;
          } else if (epsilon_n < epsilon_m){
            Y(m) = Y(m) + epsilon;
            n = m;
            m = m+1;
            // Note that nb_offsprings(n) is already set at the correct value
          } else if (epsilon_n > epsilon_m){
            Y(n) = Y(n) - epsilon;
            nb_offsprings(m) += 1;
            m = m+1;
          }
        }
      }
    }
    int idx = 0;
    for (int l = 0; l < N; l++){
      for (int j = 0; j < nb_offsprings(l); j++){
        ancestors(idx) = l;
        idx++;
      }
    }
  }
  return ancestors;
}
