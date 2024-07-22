# sample Multivariate Skewed Normal
multivar.skewed <- function(n, Pi, xi, Omega, alpha){
  # omega = vector of probabilties for mixture
  # p_mix, mu_mix, sigma_mix = parameters for 3 component mixture
  # mu, sigma = normal parameters
  # shape, rate = gamma parameters
  K <- 3
  s <- sample(1:K, size = n, replace = T, prob = Pi)
  x <- (s==1) * sn::rmsn(n=n, xi = xi[1,], Omega = Omega[,,1], alpha = alpha[1,]) +
    (s==2) * sn::rmsn(n=n, xi = xi[2,], Omega = Omega[,,2], alpha = alpha[2,]) +
    (s==3) * sn::rmsn(n=n, xi = xi[3,], Omega = Omega[,,3], alpha = alpha[3,])
  return(list(data = x, s = s))
}