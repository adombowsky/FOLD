# creates distance matrix (upper triangle only) for L multivariate normal distributions
mnorm.D <- function(mu, Sigma){
  # mu = L x p matrix of means
  # Sigma = p x p x L array of covariance matricies
  L <- nrow(mu)
  D <- matrix(0, nrow = L, ncol = L)
  for (l in 1:(L-1)){
    for (j in l:L) {
      D[l,j] = mnorm.hellinger(mu_1 = mu[l,], Sigma_1 = Sigma[,,l], mu_2 = mu[j,], Sigma_2 = Sigma[,,j])
    }
  }
  D_ut <- t(D)[lower.tri(t(D))]
  return(D_ut)
}