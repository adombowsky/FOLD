#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat rdirichlet_arma(int n, vec a){
  int K = a.size();
  vec g = zeros(K);
  mat D = zeros(n,K);
  for (int i = 0; i<n; i++){
    for (int m = 0; m < K; m++){
      vec samp = Rcpp::rgamma(1, a(m), 1);
      g(m) = as_scalar(samp);
    }
    vec dir = g/sum(g);
    D.row(i) = dir.t();
  }
  return(D);
}