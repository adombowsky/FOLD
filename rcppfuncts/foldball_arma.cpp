#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat getcGsamps(cube theta, double w, int d, int maxk){
  Environment myEnv = Environment::global_env();
  Function makecG = myEnv["makecG"];
  int S = theta.n_slices;
  int n = theta.n_rows;
  mat cG = zeros(S, n);
  for (int s = 0; s<S; s++){
    cG(s,) = makecG(theta.slice(s), w, d, maxk);
  }
  return(cG);
}
