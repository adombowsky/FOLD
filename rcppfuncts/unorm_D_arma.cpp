#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat unorm_D_arma(vec mu, vec sig){
  int n = mu.n_elem;
  mat D = zeros(n,n);
  for (int i = 0; i<(n-1); i++){
    for (int j = i+1; j<n; j++){
      D(i,j) = 1 - (sqrt( (2.0*sig(i)*sig(j))/(sig(i)*sig(i) + sig(j)*sig(j))) * exp(-(mu(i)-mu(j))*(mu(i)-mu(j))/(4.0 * (sig(i)*sig(i) + sig(j)*sig(j)))));
    }
  }
  return(D);
} 