#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat makeSigma(vec S_D, vec S_LT){
  int d = S_D.n_rows;
  mat E = zeros(d,d);
  E(trimatl_ind(size(E),-1)) = S_LT;
  mat Sigma = diagmat(S_D) + E + E.t();
  return(Sigma);
}

// [[Rcpp::export]]
cube makeSigmacube(mat Sig_diag, mat Sig_LT){
  int n = Sig_diag.n_rows;
  int d = Sig_diag.n_cols;
  cube Sigma = zeros(d,d,n);
  for (int i=0; i<n; i++){
    Sigma.slice(i) = makeSigma(Sig_diag.row(i).t(),Sig_LT.row(i).t());
  }
  return(Sigma);
}

// [[Rcpp::export]]
vec makeSigmaDet(cube Sigma, int n){
  vec dt = zeros(n);
  for (int i=0; i<n; i++){
    dt(i) = exp(0.25 * log_det_sympd(Sigma.slice(i)));
  }
  return (dt);
}
  
// [[Rcpp::export]]
mat mnorm_D_arma(mat mu, mat Sig_diag, mat Sig_LT){
  int n = mu.n_rows;
  mat D = zeros(n,n);
  cube Sigma = makeSigmacube(Sig_diag, Sig_LT);
  vec dt = makeSigmaDet(Sigma, n);
  // double dt_i = 0.0;
  // double dt_j = 0.0;
  double h_ij = 0.0;
  // hellinger distance
  for (int i = 0; i<(n-1); i++){
    //t_i = exp(0.25 * log_det_sympd(Sigma.slice(i)));
    for (int j = i+1; j<n; j++){
      //dt_j = exp(0.25 * log_det_sympd(Sigma.slice(j)));
      h_ij = Rf_fround(sqrt(1 - 
        (dt(i) * dt(j)/exp(0.5 * log_det_sympd(0.5*(Sigma.slice(i) + Sigma.slice(j)))))  * 
        exp(-(1/8.0) * dot( (mu.row(i) - mu.row(j)), inv_sympd(0.5*(Sigma.slice(i) + Sigma.slice(j))) * (mu.row(i) - mu.row(j)).t() ))),
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
vec makeHellingerAvg(cube theta, int d, int n) {
  int S = theta.n_slices;
  // initialize
  vec hgersamp = zeros(n*(n-1)/2);
  mat theta_s = zeros(n, 2*d + (d * (d-1)/2));
  mat mu = theta_s.cols(0, d-1);
  mat Sig_diag = theta_s.cols(d, 2*d-1);
  mat Sig_LT = theta_s.cols(2*d, 2*d + (d * (d-1)/2) - 1);
  mat D = zeros(n,n);
  uvec Du_ind = trimatu_ind(size(D),1);
  for (int s = 0; s<S; s++){
    theta_s = theta.slice(s);
    mu = theta_s.cols(0, d-1);
    Sig_diag = theta_s.cols(d, 2*d-1);
    Sig_LT = theta_s.cols(2*d, 2*d + (d * (d-1)/2) - 1);
    D = mnorm_D_arma(mu = mu, Sig_diag = Sig_diag, Sig_LT = Sig_LT);
    hgersamp = D(Du_ind) + hgersamp;
  }
  hgersamp =  (1/(S * 1.0)) * hgersamp;
  return(hgersamp);
}
