#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double ldmvnorm_arma(vec y, vec mu, mat Sigma, int d){
  double C = -d * 0.5 * log(2*datum::pi);
  double D = - 0.5 * log_det_sympd(Sigma);
  double E = - 0.5 * dot(y - mu, inv_sympd(Sigma) * (y-mu));
  double x = C + D + E;
  return(x);
}

// [[Rcpp::export]]
mat locmaketau(vec Pi, mat y, mat mu, mat Sigma){
  int n = y.n_rows;
  int d = y.n_cols;
  int L = mu.n_rows;
  mat tau = zeros(n,L);
  vec tau_x = zeros(L);
  for (int i = 0; i<n; i++){
    for (int l = 0; l<L; l++){
      tau_x(l) = log(Pi(l)) + ldmvnorm_arma(y.row(i).t(), mu.row(l).t(), Sigma, d);
    }
    tau.row(i) = exp(tau_x - tau_x.max()).t();
  }
  return(tau);
}
