#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double ldmvnorm_arma(arma::vec y, arma::vec mu, arma::mat Sigma, int d){
  double C = -d * 0.5 * log(2*arma::datum::pi);
  double D = - 0.5 * arma::log_det_sympd(Sigma);
  double E = - 0.5 * arma::dot(y - mu, inv_sympd(Sigma) * (y-mu));
  double x = C + D + E;
  return(x);
}

// [[Rcpp::export]]
double ldunorm_arma(double y, double mu, double sigma_sq){
  double C = -0.5 * log(2*arma::datum::pi);
  double D = - 0.5 * log(sigma_sq);
  double E = - 0.5 * pow(y-mu,2.0)/sigma_sq;
  double x = C + D + E;
  return(x);
}

// [[Rcpp::export]]
arma::mat maketau(arma::vec Pi, arma::mat y, arma::mat mu, arma::cube Sigma){
  int n = y.n_rows;
  int d = y.n_cols;
  int L = mu.n_rows;
  arma::mat tau = arma::zeros(n,L);
  arma::vec tau_x = arma::zeros(L);
  for (int i = 0; i<n; i++){
    for (int l = 0; l<L; l++){
      tau_x(l) = log(Pi(l)) + ldmvnorm_arma(y.row(i).t(), mu.row(l).t(), Sigma.slice(l), d);
    }
    tau.row(i) = arma::exp(tau_x - tau_x.max()).t();
  }
  return(tau);
}

// [[Rcpp::export]]
arma::mat makeutau(arma::vec Pi, arma::vec y, arma::vec mu, arma::vec Sigma){
  int n = y.n_rows;
  int L = mu.n_rows;
  arma::mat tau = arma::zeros(n,L);
  arma::vec tau_x = arma::zeros(L);
  for (int i = 0; i<n; i++){
    for (int l = 0; l<L; l++){
      tau_x(l) = log(Pi(l)) + ldunorm_arma(y(i), mu(l), Sigma(l));
    }
    tau.row(i) = arma::exp(tau_x - tau_x.max()).t();
  }
  return(tau);
}

// [[Rcpp::export]]
arma::mat lmaketau(arma::vec Pi, arma::mat y, arma::mat mu, arma::mat Sigma){
  int n = y.n_rows;
  int d = y.n_cols;
  int L = mu.n_rows;
  arma::mat tau = arma::zeros(n,L);
  arma::vec tau_x = arma::zeros(L);
  for (int i = 0; i<n; i++){
    for (int l = 0; l<L; l++){
      tau_x(l) = log(Pi(l)) + ldmvnorm_arma(y.row(i).t(), mu.row(l).t(), Sigma, d);
    }
    tau.row(i) = arma::exp(tau_x - tau_x.max()).t();
  }
  return(tau);
}
