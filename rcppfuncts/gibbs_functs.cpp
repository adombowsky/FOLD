#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
using namespace Rcpp;
arma::mat cont_sij_and_yij(int n, int K, arma::vec y_j, arma::mat theta_j,
                           arma::ivec s_j, arma::uvec na_j, arma::mat l_eta){
  arma::vec pr = arma::zeros(K);
  for (int i = 0; i<n; i++){
   if (na_j(i)) {
     y_j(i) = theta_j(s_j(i), 0) + sqrt(theta_j(s_j(i), 2)) * arma::randn() 
     }
   else {
     y_j(i) = y_j(i)
     }
    for (int h=0; h<K; h++){
      pr(h) = l_eta(i,h) + 
        dnorm(x = y_j(i), mean = theta_j(0, h), sd = sqrt(theta_j(1,h)), 
                    log = TRUE) 
      }
    arma::vec prob = exp(pr - max(pr));
    s_j(i) = sample(x = 1:K,  size = 1, prob = prob);
  }
  arma::mat Z = cbind(y_j, s_j);
  return(Z);
}


p <- rep(0, K)
  for (h in 1:K){
    if (type == "c"){ # continuous
      p[h] <- log(eta_i[h]) + 
        dnorm(x = y_ij, mean = theta_j[1,h], sd = sqrt(theta_j[2,h]), 
              log = T)
    }
    else if (type == "b"){ # binary
      p[h] <- log(eta_i[h]) + 
        dbinom(x = y_ij, size = 1, prob = theta_j[h], log = T)
    }
  }
  p <- exp(p - max(p))
    s_ij <- sample(x = 1:K, size = 1, replace = T, prob = p)


for (i in 1:n){
  if (na_matrix[i, j]) {
    y_ij <- sample.y_ij(theta_j = theta[[j]][r, , ], 
                        type = "c", s_ij = s[i,j])
  } else{
    y_ij <- y[i, j] 
  }
  s[i, j] <- sample.s_ij(y_ij = y_ij, type = "c", eta_i = eta[r-1, , i], 
                         theta_j = t(theta[[j]][r, , ]))
}