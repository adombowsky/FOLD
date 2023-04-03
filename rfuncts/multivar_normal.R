# sample Multivariate Skewed Normal
multivar.normal <- function(n, Pi, xi, Omega){
  # omega = vector of probabilties for mixture
  # p_mix, mu_mix, sigma_mix = parameters for 3 component mixture
  # mu, sigma = normal parameters
  # shape, rate = gamma parameters
  K <- 3
  s <- sample(1:K, size = n, replace = T, prob = Pi)
  x <- (s==1) * mvtnorm::rmvnorm(n=n, mean = xi[1,], sigma = Omega[,,1]) +
    (s==2) * mvtnorm::rmvnorm(n=n, mean = xi[2,], sigma = Omega[,,2]) +
    (s==3) * mvtnorm::rmvnorm(n=n, mean = xi[3,], sigma = Omega[,,3])
  return(list(data = x, s = s))
}