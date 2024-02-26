#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double risk_cpp(arma::vec c, arma::mat Delta, double omega) {
  int n = c.n_rows;
  double risk = 0.0;
  for (int i = 0; i<(n-1); i++){
    for (int j = i+1; j<n; j++){
      risk = risk + (c(i)==c(j))*Delta(i,j) + (c(i) != c(j))*omega*(1-Delta(i,j));
    }
    }
  return(risk);
}

// [[Rcpp::export]]
double rand_index(arma::rowvec c1, arma::rowvec c2) {
  int n = c1.n_cols;
  arma::mat RI_mat = arma::zeros(n,n);
  for (int i=0; i<(n-1); i++) {
    for (int j=i+1; j<n; j++) {
      RI_mat(i,j) = ((c1(i) == c1(j)) & (c2(i) == c2(j))) + ((c1(i) != c1(j)) & (c2(i) != c2(j)));
    }
  }
  return(arma::accu(RI_mat)/(0.5 * n * (n-1)));
}

// [[Rcpp::export]]
arma::vec risk_matrix_cpp(arma::mat c, arma::mat Delta, double omega) {
  int M = c.n_rows;
  arma::vec risk = arma::zeros(M,1);
  for (int m = 0; m < M; m++) {
    risk(m) = risk_cpp(c.row(m).t(), Delta, omega);
  }
  return(risk);
}