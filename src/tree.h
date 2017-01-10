#ifndef _INCL_BushT_
#define _INCL_BushT_
#include <RcppEigen.h>
using namespace Rcpp;

class Tree
{
  public:
    Tree(int N, int M, int dimx);
  virtual ~Tree();
  // methods
  void reset();
  void init(NumericMatrix x);
  void insert(NumericMatrix x, IntegerVector a);
  void prune(IntegerVector o);
  void update(NumericMatrix x, IntegerVector a);
  NumericMatrix retrieve_xgeneration(int lag);
  void double_size();
  NumericMatrix get_path(int n);
  // attributes
  // number of particles
  int N;
  // upper bound on total size of the tree
  int M;
  // dimension of each x-particle
  int dimx;
  // number of time "insert" has been called
  int nsteps;
  // vector of ancestor indices
  IntegerVector a_star;
  // vector of offspring counts
  IntegerVector o_star;
  // vector of particles
  NumericMatrix x_star;
  // vector of leaf indices
  IntegerVector l_star;
  //

};

#endif
