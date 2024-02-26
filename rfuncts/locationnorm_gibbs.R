# other functions
ldMvn <- function(X,mu,Sigma,p) {
  rooti <- backsolve(chol(Sigma),diag(p))
  quads <- colSums((crossprod(rooti,(X-mu)))^2)
  return(-(p/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads)
}



loc_mvnorm_gibbs <- function(S, y, L, alpha=rep(1/L,L), mu_0, Sigma_0, Sigma, stops = 1){
  # the DEFINITIVE gibbs sampler for a multivariate hierarchical Bayesian GMM to compute Hellinger
  
  # packages
  require(LaplacesDemon)
  require(mclust)
  
  # Rcpp functions
  Rcpp::sourceCpp("rcppfuncts/loc_clusterlabs_arma.cpp")
  
  
  ## S = number of iterations, y = data
  ## L = number of kernels
  ## alpha = Dirichlet prior parameter
  ## w = mean location parameter
  ## kappa = mean scale parameter
  ## r = inverse wishart degrees of freedom, r > p-1
  ## C = inverse wishart scale matrix, C > 0
  
  # prelims
  n <- nrow(y)
  p <- ncol(y)
  # initialization
  Pi <- alpha/sum(alpha) # weights
  mu <- matrix(0, nrow = L, ncol = p) # component means
  #Sigma <- array(diag(1,p), dim = c(p,p,L))
  #Sigma_diag <- matrix(1, nrow = L, ncol = p) # component variances
  #Sigma_LT <- matrix(1, nrow = L, ncol = choose(p,2)) # lower triangles
  Sigma_0_inv <- solve(Sigma_0)
  Sigma_inv <- solve(Sigma)
  z <- matrix(sample(1:L, size = n*S, replace = T), nrow = S, ncol = n ) # cluster labelså
  tau <- rep(0,L) # for sampling cluster labels
  theta <- list() # for Hellinger distance
  E <- matrix(0, nrow = p, ncol = p)
  # configuring stops
  stops <- (1:(S/stops)) * stops
  # sampling
  print("Sampling")
  for (s in 2:S){
    # compute cluster sizes
    cluster.sizes <- as.vector(table(factor(z[s-1, ], levels = 1:L)))
    # sample pi
    alpha_star <- alpha + cluster.sizes
    Pi <- LaplacesDemon::rdirichlet(n=1, alpha=alpha_star)
    # sample mu and Sigma
    for (l in 1:L) {
      if (cluster.sizes[l] == 0) { # for empty clusters, sample from the prior
        mu[l, ] <- mu_0 + chol(Sigma_0) %*% rnorm(p)
      } 
      else if (cluster.sizes[l] == 1) {
        ## sampling mean and covariance
        y_l = y[z[s-1,] == l, ]
        n_l = sum(z[s-1,]==l)
        y_l_bar = y_l
        Sigma_l <- solve(Sigma_0_inv + n_l * Sigma_inv)
        mu_l <- Sigma_l %*% (Sigma_0_inv %*% mu_0 + n_l * Sigma_inv %*% y_l_bar)
        ## mean sampled
        mu[l,] <- mu_l + chol(Sigma_l) %*% rnorm(p)
      }
      else{
        ## sampling mean and variance
        y_l = y[z[s-1,] == l, ]
        n_l = sum(z[s-1,]==l)
        y_l_bar = colMeans(y_l)
        Sigma_l <- solve(Sigma_0_inv + n_l * Sigma_inv)
        mu_l <- Sigma_l %*% (Sigma_0_inv %*% mu_0 + n_l * Sigma_inv %*% y_l_bar)
        ## mean sampled
        mu[l,] <- mu_l + chol(Sigma_l) %*% rnorm(p)
      }
    }
    # sample cluster labels
    tau <- locmaketau(Pi = Pi, y = y, mu = mu, Sigma = Sigma)
    z[s,] <- apply(tau, 1, function(x) sample(L, size = 1, replace = T, prob = x))
    # sub-sampling for theta samples
    z_mat <- model.matrix(~ -1 + cl, data = data.frame(cl = factor(z[s,], levels = 1:L)))
    total_pars <- cbind(mu)
    theta[[s]] <- z_mat %*% total_pars
    
    # print stops
    if (s %in% stops){
      print(s)
    }
  }
  return(list(theta = theta, z=z))
}