#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat locnorm_D_arma(mat mu, mat Sigma){
  int n = mu.n_rows;
  mat D = zeros(n,n);
  mat Sigma_inv = inv_sympd(Sigma);
  // double dt_i = 0.0;
  // double dt_j = 0.0;
  double h_ij = 0.0;
  // hellinger distance
  for (int i = 0; i<(n-1); i++){
    //t_i = exp(0.25 * log_det_sympd(Sigma.slice(i)));
    for (int j = i+1; j<n; j++){
      //dt_j = exp(0.25 * log_det_sympd(Sigma.slice(j)));
      h_ij = Rf_fround(sqrt(1 - 
        exp(-(1/8.0) * dot( (mu.row(i) - mu.row(j)), Sigma_inv * (mu.row(i) - mu.row(j)).t() ))),
        5);
      if (h_ij >= 0){
        D(i,j) = h_ij;
      } else {
        D(i,j) = 0.0;
      }
    }
  }
  return(D);
} 

// [[Rcpp::export]]
vec locHellingerAvg(cube theta, mat Sigma, int d, int n) {
  int S = theta.n_slices;
  // initialize
  vec hgersamp = zeros(n*(n-1)/2);
  mat mu = zeros(n, d);
  mat D = zeros(n,n);
  uvec Du_ind = trimatu_ind(size(D),1);
  for (int s = 0; s<S; s++){
    mu = theta.slice(s);
    D = locnorm_D_arma(mu = mu, Sigma = Sigma);
    hgersamp = D(Du_ind) + hgersamp;
  }
  hgersamp =  (1/(S * 1.0)) * hgersamp;
  return(hgersamp);
}
