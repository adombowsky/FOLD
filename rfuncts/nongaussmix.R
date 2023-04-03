# sample Multivariate Skewed Normal
multivar.gmm <- function(n, p, Mu, Sigma){
  s <- sample(1:2, size = n, replace = T, prob = c(p,1-p))
  x <- (s==1) * rmvnorm(n = n, mean = Mu[1,], sigma = Sigma[,,1]) +
    (s==2) * rmvnorm(n=n, mean = Mu[2,], sigma = Sigma[,,2])
  return(x)
}
nongauss <- function(n, Pi, xi, Omega, alpha, Mu_1, Sigma_1, Mu, Sigma, p){
  # omega = vector of probabilties for mixture
  # p_mix, mu_mix, sigma_mix = parameters for 3 component mixture
  # mu, sigma = normal parameters
  # shape, rate = gamma parameters
  require(sn)
  K <- 4
  s <- sample(1:K, size = n, replace = T, prob = Pi)
  x <- (s==1) * rmsn(n=n, xi = xi[1,], Omega = Omega[,,1], alpha = alpha[1,]) +
    (s==2) * rmsn(n=n, xi = xi[2,], Omega = Omega[,,2], alpha = alpha[2,]) +
    (s==3) * rmvnorm(n=n, mean = Mu_1, sigma = Sigma_1) +
    (s==4) * multivar.gmm(n = n, p = p, Mu = Mu, Sigma = Sigma)
  return(list(data = x, s = s))
}