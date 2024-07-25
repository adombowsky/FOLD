#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat makeSigma(arma::vec S_D, arma::vec S_LT){
  int d = S_D.n_rows;
  arma::mat E = arma::zeros(d,d);
  E(trimatl_ind(size(E),-1)) = S_LT;
  arma::mat Sigma = arma::diagmat(S_D) + E + E.t();
  return(Sigma);
}


arma::cube makeSigmacube(arma::mat Sig_diag, arma::mat Sig_LT){
  int n = Sig_diag.n_rows;
  int d = Sig_diag.n_cols;
  arma::cube Sigma = arma::zeros(d,d,n);
  for (int i=0; i<n; i++){
    Sigma.slice(i) = makeSigma(Sig_diag.row(i).t(),Sig_LT.row(i).t());
  }
  return(Sigma);
}


arma::vec makeSigmaDet(arma::cube Sigma, int n){
  arma::vec dt = arma::zeros(n);
  for (int i=0; i<n; i++){
    dt(i) = exp(0.25 * arma::log_det_sympd(Sigma.slice(i)));
  }
  return (dt);
}

// [[Rcpp::export]]
arma::mat mnorm_D_arma(arma::mat mu, arma::mat Sig_diag, arma::mat Sig_LT){
  int n = mu.n_rows;
  arma::mat D = arma::zeros(n,n);
  arma::cube Sigma = makeSigmacube(Sig_diag, Sig_LT);
  arma::vec dt = makeSigmaDet(Sigma, n);
  double h_ij = 0.0;
  // hellinger distance
  for (int i = 0; i<(n-1); i++){
    for (int j = i+1; j<n; j++){
      h_ij = Rf_fround(sqrt(1 -
        (dt(i) * dt(j)/exp(0.5 * arma::log_det_sympd(0.5*(Sigma.slice(i) + Sigma.slice(j)))))  *
        exp(-(1/8.0) * dot( (mu.row(i) - mu.row(j)), arma::inv_sympd(0.5*(Sigma.slice(i) + Sigma.slice(j))) * (mu.row(i) - mu.row(j)).t() ))),
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
arma::mat unorm_D_arma(arma::vec mu, arma::vec sigma){
  int n = mu.n_rows;
  arma::mat D = arma::zeros(n,n);
  double h_ij = 0.0;
  // hellinger distance
  for (int i = 0; i<(n-1); i++){
    for (int j = i+1; j<n; j++){
      h_ij = Rf_fround(sqrt(1 -
      pow(2.0 * sigma(i) * sigma(j)/(pow(sigma(i),2) + pow(sigma(j),2)), 0.5) *
      exp(-0.25 * pow(mu(i)-mu(j), 2.0)/(pow(sigma(i),2) + pow(sigma(j),2)))),
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
arma::vec makeHellingerAvg(arma::cube theta, int d, int n) {
  int S = theta.n_slices;
  // initialize
  arma::vec hgersamp = arma::zeros(n*(n-1)/2);
  arma::mat theta_s = arma::zeros(n, 2*d + (d * (d-1)/2));
  arma::mat mu = theta_s.cols(0, d-1);
  arma::mat Sig_diag = theta_s.cols(d, 2*d-1);
  arma::mat Sig_LT = theta_s.cols(2*d, 2*d + (d * (d-1)/2) - 1);
  arma::mat D = arma::zeros(n,n);
  arma::uvec Du_ind = arma::trimatu_ind(size(D),1);
  for (int s = 0; s<S; s++){
    theta_s = theta.slice(s);
    mu = theta_s.cols(0, d-1);
    Sig_diag = theta_s.cols(d, 2*d-1);
    Sig_LT = theta_s.cols(2*d, 2*d + (d * (d-1)/2) - 1);
    D = mnorm_D_arma(mu, Sig_diag, Sig_LT);
    hgersamp = D(Du_ind) + hgersamp;
  }
  hgersamp =  (1/(S * 1.0)) * hgersamp;
  return(hgersamp);
}

// [[Rcpp::export]]
arma::vec makeuHellingerAvg(arma::cube theta, int n) {
  int S = theta.n_slices;
  // initialize
  arma::vec hgersamp = arma::zeros(n*(n-1)/2);
  arma::mat theta_s = arma::zeros(n, 2);
  arma::vec mu = theta_s.col(0);
  arma::vec sigma = arma::sqrt(theta_s.col(1));
  arma::mat D = arma::zeros(n,n);
  arma::uvec Du_ind = arma::trimatu_ind(size(D),1);
  for (int s = 0; s<S; s++){
    theta_s = theta.slice(s);
    mu = theta_s.col(0);
    sigma = arma::sqrt(theta_s.col(1));
    D = unorm_D_arma(mu, sigma);
    hgersamp = D(Du_ind) + hgersamp;
  }
  hgersamp =  (1/(S * 1.0)) * hgersamp;
  return(hgersamp);
}
