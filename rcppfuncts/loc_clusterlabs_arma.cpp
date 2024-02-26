#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double ldmvnorm_arma(arma::vec y, arma::vec mu, arma::mat Sigma, int d){
  double C = -d * 0.5 * log(2*arma::datum::pi);
  double D = - 0.5 * arma::log_det_sympd(Sigma);
  double E = - 0.5 * dot(y - mu, inv_sympd(Sigma) * (y-mu));
  double x = C + D + E;
  return(x);
}

// [[Rcpp::export]]
int sample_arma(arma::vec probs) {
  int K = probs.n_rows;
  IntegerVector clusts = Range(1,K);
  IntegerVector samp = RcppArmadillo::sample(clusts, 1, TRUE, probs);
  int s = samp(0);
  return(s);
}

// [[Rcpp::export]]
arma::mat locmaketau(arma::vec Pi, arma::mat y, arma::mat mu, arma::mat Sigma){
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

// [[Rcpp::export]]
arma::vec gibbs_C_location(arma::vec c,
                          int L,
                          arma::mat mu,
                          arma::mat X,
                          arma::vec alpha,
                          double sigma_sq){
  int n = X.n_rows;
  int d = X.n_cols;
  int n_k = n;
  arma::vec lp = arma::zeros(L,1);
  arma::vec p = arma::zeros(L,1);
  for (int i=0; i<n; i++) {
    for (int k=0; k<L; k++) {
      if (c(i)==k+1) {
        n_k = arma::sum(c==k+1)-1;
      } else {
        n_k = arma::sum(c==k+1);
      }
      lp(k) = log(n_k + (alpha(k)/(L*1.0))) + ldmvnorm_arma(X.row(i).t(), mu.row(k).t(), sigma_sq * arma::eye(d,d), d);
    }
    p = arma::exp(lp - arma::max(lp));
    c(i) = sample_arma(p/sum(p));
  }
  return(c);
}

