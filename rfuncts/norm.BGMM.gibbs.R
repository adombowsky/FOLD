
sample.labels <- function(x, omega, mu, sigma_sq, L){
  prbs <- log(omega) + dnorm(x, mean = mu, sd = sqrt(sigma_sq), log = T)
  prbs <- exp(prbs - max(prbs))
  Z <- sample(1:L, size = 1, replace = T, prob = prbs)
  return(Z)
}

sample.mu.sigma <- function(x, y, Z){
  # get parameters
  mu_0 <- x[1]
  kappa <- x[2]
  a_sq <- x[3]
  b_sq <- x[4]
  l <- x[5]
  y_l <- y[Z == l]
  n_l <- length(y_l)
  # fit
  kappa_hat <- 1/((1/kappa) + n_l)
  mu_hat <- kappa_hat * ( mu_0/kappa + n_l * ifelse(n_l > 0, mean(y_l), 0))
  a_hat <- a_sq + n_l/2
  b_hat <- b_sq + 0.5 * (ifelse(n_l > 0, sum( (y_l - mean(y_l))^2), 0) + 
                              ifelse(n_l > 0, (n_l/(1 + kappa * n_l)) * (mean(y_l) - mu_0)^2, 0 ))
  sigma_sq <- 1/rgamma(n = 1, shape = a_hat, rate = b_hat)
  mu <- rnorm(n = 1, mean = mu_hat, sd = sqrt(kappa_hat * sigma_sq))
  return(c(mu, sigma_sq))
}

ugmm.overfitted.gibbs <- function(R, L, y, priors, alpha, stops) {
  # at long last, a DEFINITIVE function for fitting the univariate overfitted GMM
  # note that lots of components will be empty, emphasis should be on speed
  require(Rcpp)
  require(RcppArmadillo)
  
  
  print("Compiling C++ Functions")
  sourceCpp("rcppfuncts/rdirichlet_arma.cpp")
  
  n <- length(y)
  Z <- matrix(sample(1:L, size = n * R, replace = T, prob = rep(1/L, L)), nrow = R, ncol = n)
  mu <- matrix(0, nrow = R, ncol = L)
  sigma_sq <- matrix(1, nrow = R, ncol = L)
  omega <- matrix(1/L, nrow = R, ncol = L)
  
  # configuring stops
  stops <- (1:(R/stops)) * stops
  for (r in 2:R){
    # update labels
    Z[r,] <- sapply(y, sample.labels, omega = omega[r-1,], mu = mu[r-1,], sigma_sq = sigma_sq[r-1,], L = L)
    # update components
    pars <- apply(priors, 1, sample.mu.sigma, y = y, Z = Z[r,])
    mu[r, ] <- pars[1,]
    sigma_sq[r, ] <- pars[2,]
    # update weights
    n_l <- table(factor(Z[r,], levels = 1:L))
    omega[r, ] <- rdirichlet_arma(n = 1, a = alpha + n_l)
  # print stops
    if (r %in% stops){
      print(r)
    } 
    }
  # returning
  return(list(Z = Z, mu = mu, sigma_sq = sigma_sq, omega = omega))
  
  
}