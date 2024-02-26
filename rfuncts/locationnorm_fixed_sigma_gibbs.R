loc_fsig_mvnorm_gibbs_cpp <- function(S, y, L, alpha=rep(1/L,L), mu_0, Sigma_0, sigma_sq, stops = 1){
  # the DEFINITIVE gibbs sampler for a multivariate hierarchical Bayesian GMM to compute Hellinger
  
  # packages
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
  mu <- matrix(rnorm(L*p), nrow = L, ncol = p) # component means
  #Sigma <- array(diag(1,p), dim = c(p,p,L))
  #Sigma_diag <- matrix(1, nrow = L, ncol = p) # component variances
  #Sigma_LT <- matrix(1, nrow = L, ncol = choose(p,2)) # lower triangles
  Sigma_0_inv <- solve(Sigma_0)
  Sigma <- sigma_sq * diag(p)
  Sigma_inv <- solve(Sigma)
  z <- matrix(sample(1:L, size = n*S, replace = T), nrow = S, ncol = n ) # cluster labels
  theta <- list() # for Hellinger distance
  E <- matrix(0, nrow = p, ncol = p)
  # configuring stops
  stops <- (1:(S/stops)) * stops
  # sampling
  print("Sampling")
  for (s in 2:S){
    # sample mu
    cluster.sizes <- as.vector(table(factor(z[s-1, ], levels = 1:L)))
    for (l in 1:L) {
      if (cluster.sizes[l] == 0) { # for empty clusters, sample from the prior
        #Sigma[,,l] <- LaplacesDemon::rinvwishart(nu = r, S = C)
        #Sigma_diag[l, ] <- diag(Sigma[,,l])
        #Sigma_LT[l, ] <- Sigma[,,l][lower.tri(Sigma[,,l], diag = F)]
        mu[l, ] <- mu_0 + chol(Sigma_0) %*% rnorm(p)
      } 
      else if (cluster.sizes[l] == 1) {
        ## sampling mean and covariance
        y_l = y[z[s-1,] == l, ]
        n_l = sum(z[s-1,]==l)
        y_l_bar = y_l
        Sigma_l <- solve(Sigma_0_inv + n_l * Sigma_inv)
        mu_l <- Sigma_l %*% (Sigma_0_inv %*% mu_0 + n_l * Sigma_inv %*% y_l_bar)
        #C_star <- C +  (n_l * kappa/(n_l + kappa)) * tcrossprod(y_l_bar - w)
        #Sigma[,,l] <- rinvwishart(nu = n_l + r, S = C_star)
        #Sigma_diag[l, ] <- diag(Sigma[,,l])
        #Sigma_LT[l, ] <- Sigma[,,l][lower.tri(Sigma[,,l], diag = F)]
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
        #C_star <- C +  (n_l * kappa/(n_l + kappa)) * tcrossprod(y_l_bar - w)
        #Sigma[,,l] <- rinvwishart(nu = n_l + r, S = C_star)
        #Sigma_diag[l, ] <- diag(Sigma[,,l])
        #Sigma_LT[l, ] <- Sigma[,,l][lower.tri(Sigma[,,l], diag = F)]
        ## mean sampled
        mu[l,] <- mu_l + chol(Sigma_l) %*% rnorm(p)
      }
    }
    # sample cluster labels
    z[s,] <- gibbs_C_location(c=matrix(z[s-1,]),
                              L=L,
                              mu=mu,
                              X=y,
                              alpha=alpha,
                              sigma_sq = sigma_sq) 
  
    theta[[s]] <- mu[z[s,],]
    
    # print stops
    if (s %in% stops){
      print(s)
    }
  }
  return(list(theta = theta, z=z))
}