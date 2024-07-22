#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double dmv(rowvec X, rowvec mu, mat Sigma){
  int p = X.n_cols;
  double K = pow(2.0 * M_PI, - p/2.0);
  double g = K * exp(-0.5 * log_det_sympd(Sigma)- 0.5 * dot(X - mu, inv_sympd(Sigma) * (X - mu).t()));
  return(g);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double hgauss(rowvec mu_1, rowvec mu_2, mat Sigma_1, mat Sigma_2){
  double h = 
    Rf_fround(sqrt(1 - exp(0.25 * log_det_sympd(Sigma_1) + 
    0.25 * log_det_sympd(Sigma_2) -
    0.5 * log_det_sympd( (Sigma_1 + Sigma_2)/2.0) - 
    (1/8.0) * dot((mu_1 - mu_2), inv_sympd(0.5 * (Sigma_1 + Sigma_2)) * (mu_1 - mu_2).t()))), 5);
  return(h);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat makeh(mat mu, cube Sigma){
  int L = mu.n_rows;
  mat h = zeros(L,L);
  for (int l = 0; l<(L-1); l++){
    for (int m = l+1; m<L; m++){
      h(l,m) = hgauss(mu.row(l), mu.row(m), Sigma.slice(l), Sigma.slice(m));
    }
  }
  return(h);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double Deltapair(rowvec X_i, rowvec X_j, mat mu, cube Sigma, vec pi, mat h){
  int L = mu.n_rows;
  mat q = zeros(L, L);
  double f_i = pi(0)*dmv(X_i, mu.row(0), Sigma.slice(0))+
    pi(1)*dmv(X_i, mu.row(1), Sigma.slice(1))+
    pi(2)*dmv(X_i, mu.row(2), Sigma.slice(2));
  double f_j = pi(0)*dmv(X_j, mu.row(0), Sigma.slice(0))+
    pi(1)*dmv(X_j, mu.row(1), Sigma.slice(1))+
    pi(2)*dmv(X_j, mu.row(2), Sigma.slice(2));
  for (int l = 0; l<(L-1); l++){
    for (int m = l+1; m<L; m++){
      q(l,m) = h(l,m) *
        pi(l) * pi(m) * (dmv(X_i, mu.row(l), Sigma.slice(l)) * dmv(X_j, mu.row(m), Sigma.slice(m)) +
        dmv(X_i, mu.row(m), Sigma.slice(m)) * dmv(X_j, mu.row(l), Sigma.slice(l)))/(f_i * f_j);
    }
  }
  double D_ij = accu(q);
  return(D_ij);
}

//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat makeOracle(mat X, mat mu, cube Sigma, vec pi, mat h){
  int n = X.n_rows;
  mat Deltastar = zeros(n,n);
  for (int i = 0; i<(n-1); i++){
    for (int j = i+1; j<n; j++){
      Deltastar(i,j) = Deltapair(X.row(i), X.row(j), mu, Sigma, pi, h);
    }
  }
  return(Deltastar);
}



