#include <RcppEigen.h>
#include "tree.h"
using namespace Rcpp;
using namespace std;

Tree::Tree(int N, int M, int dimx) : N(N), M(M), dimx(dimx) {
  if (this->M < 3*this->N){
    this->M = 3*this->N;
  }
  this->nsteps = 0;
  a_star = IntegerVector(this->M);
  o_star = IntegerVector(this->M);
  x_star = NumericMatrix(this->dimx, this->M);
  l_star = IntegerVector(this->N);
  this->reset();
}

Tree::~Tree(){
  //
}

void Tree::reset(){
  this->nsteps = 0;
  std::fill(a_star.begin(), a_star.end(), 0);
  std::fill(o_star.begin(), o_star.end(), 0);
}

void Tree::init(NumericMatrix x_0){
  for (int i = 0; i < N; i ++){
    a_star(i) = -2;
    x_star(_,i) = x_0(_,i);
    l_star(i) = i;
  }
}

void Tree::insert(NumericMatrix x, IntegerVector a){
  nsteps ++;
  // b_t <- gather(l_star, a) // indices of the parents of new generation
  IntegerVector b(N);
  for (int i = 0; i < N; i ++){
    b(i) = l_star(a(i));
  }
  // z_star <- transform prefix sum(o_star, 1_0)
  // l_star <- lower bound (z_star, [1,...,N])
  int slot = 0;
  int i = 0;
  while (slot < N && i < M){
    if (o_star(i) == 0){
      l_star(slot) = i;
      slot ++;
    }
    i++;
  } // here we can test whether slot == N, ie whether we found enough slots
  if (slot < N){
    // if not double the size of vectors
    int old_M = this->M;
    double_size();
    for (i = slot; i < N; i++){
      l_star(i) = old_M + i - slot;
    }
  }
  // a_star <- scatter(b_t, l_star)
  // x_star <- scatter(x_t, l_star)
  for (int i = 0; i < N; i ++){
    a_star(l_star(i)) = b(i);
    x_star(_,l_star(i)) = x(_,i);
  }
}

void Tree::prune(IntegerVector o){
  // o_star <- scatter(o_t, l_star)
  for (int i = 0; i < N; i ++){
    o_star(l_star(i)) = o(i);
  }
  //
    int j, new_j;
  for (int i = 0; i < N; i ++){
    j = l_star(i);
    while (j >= 0 && o_star(j) == 0){
      j = a_star(j);
      if (j >= 0){
        o_star(j) = o_star(j) - 1;
      }
    }
  }
}

void Tree::update(NumericMatrix x, IntegerVector a){
  // convert ancestor vector to offspring vector
  IntegerVector o(N);
  std::fill(o.begin(), o.end(), 0);
  for (int i = 0; i < N; i++){
    o(a(i)) ++;
  }
  // prune tree
  this->prune(o);
  // insert new generation
  this->insert(x, a);
}

void Tree::double_size(){
  int old_M = this->M;
  this->M = 2*old_M;
  //  cout << "doubling size of tree to " << this->M << endl;
  IntegerVector new_a_star(this->M);
  IntegerVector new_o_star(this->M);
  NumericMatrix new_x_star(this->dimx, this->M);
  for (int i = 0; i < old_M; i++){
    new_a_star(i) = a_star(i);
    new_x_star(_,i) = x_star(_,i);
    new_o_star(i) = o_star(i);
  }
  this->a_star = new_a_star;
  this->x_star = new_x_star;
  this->o_star = new_o_star;
}

NumericMatrix Tree::get_path(int n){
  NumericMatrix path(dimx, this->nsteps + 1);
  int j = this->l_star(n);
  path(_,this->nsteps) = this->x_star(_,j);
  int step = this->nsteps - 1;
  while (j >= 0 && step >= 0){
    j = this->a_star(j);
    path(_,step) = this->x_star(_,j);
    step --;
  }
  return path;
}

NumericMatrix Tree::retrieve_xgeneration(int lag){
  NumericMatrix xgeneration(dimx, N);
  for (int i_particle = 0; i_particle < N; i_particle ++){
    int j = this->l_star(i_particle);
    for (int i_step = 0; i_step < lag; i_step++){
      j = this->a_star(j);
    }
    xgeneration(_,i_particle) = this->x_star(_,j);
  }
  return xgeneration;
}

// Added a new method to reconstruct tree from a list of attributes:
// Necessary to rebuild tree after loading attributes from RDS file of a partial run save
void Tree::reconstruct(int N, int M, int d, int n, IntegerVector a, IntegerVector o, NumericMatrix x, IntegerVector l){
  this->N = N;
  this->M = M;
  this->dimx = d;
  this->nsteps = n;
  this->a_star = a;
  this->o_star = o;
  this->x_star = x;
  this->l_star = l;
}

RCPP_MODULE(module_tree) {
  class_<Tree>( "Tree" )
  .constructor<int,int,int>()
  .field( "N", &Tree::N)
  .field( "M", &Tree::M)
  .field( "dimx", &Tree::dimx)
  .field( "nsteps", &Tree::nsteps)
  .field( "a_star", &Tree::a_star)
  .field( "o_star", &Tree::o_star)
  .field( "x_star", &Tree::x_star)
  .field( "l_star", &Tree::l_star)
  .method( "init", &Tree::init)
  .method( "insert", &Tree::insert)
  .method( "prune", &Tree::prune)
  .method( "update", &Tree::update)
  .method( "get_path", &Tree::get_path)
  .method( "retrieve_xgeneration", &Tree::retrieve_xgeneration)
  .method( "reconstruct", &Tree::reconstruct)
  ;
}
