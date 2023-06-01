# other functions
ldMvn <- function(X,mu,Sigma,p) {
  rooti <- backsolve(chol(Sigma),diag(p))
  quads <- colSums((crossprod(rooti,(X-mu)))^2)
  return(-(p/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads)
}



mvnorm_gibbs <- function(S, y, L, alpha=rep(1/L,L), w, kappa, r, C, stops = 1){
  # the DEFINITIVE gibbs sampler for a multivariate hierarchical Bayesian GMM to compute Hellinger
  # from MacLachlan and Peel, pg. 123
  
  # packages
  require(LaplacesDemon)
  require(mclust)
  
  # Rcpp functions
  Rcpp::sourceCpp("rcppfuncts/clusterlabs_arma.cpp")
  
  # r functions
  #source("rfuncts/mnorm_hellinger.R")
  #source("rfuncts/mnorm_D.R")
  
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
  Sigma <- array(diag(1,p), dim = c(p,p,L))
  Sigma_diag <- matrix(1, nrow = L, ncol = p) # component variances
  Sigma_LT <- matrix(1, nrow = L, ncol = choose(p,2)) # lower triangles
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
        Sigma[,,l] <- LaplacesDemon::rinvwishart(nu = r, S = C)
        Sigma_diag[l, ] <- diag(Sigma[,,l])
        Sigma_LT[l, ] <- Sigma[,,l][lower.tri(Sigma[,,l], diag = F)]
        mu[l, ] <- w + (chol(Sigma[,,l])/sqrt(kappa)) %*% rnorm(p)
      } 
      else if (cluster.sizes[l] == 1) {
        ## sampling mean and covariance
        y_l = y[z[s-1,] == l, ]
        n_l = sum(z[s-1,]==l)
        y_l_bar = y_l
        ## covariance matrix sampled
        # V <- y_l - matrix(rep(y_l_bar,n_l), nrow = n_l, byrow = T)
        # V <- crossprod(data.matrix(V))
        # C_star <- C +  V + (n_l * r/(n_l + r)) * tcrossprod(y_l_bar - w)
        C_star <- C +  (n_l * kappa/(n_l + kappa)) * tcrossprod(y_l_bar - w)
        Sigma[,,l] <- rinvwishart(nu = n_l + r, S = C_star)
        Sigma_diag[l, ] <- diag(Sigma[,,l])
        Sigma_LT[l, ] <- Sigma[,,l][lower.tri(Sigma[,,l], diag = F)]
        ## mean sampled
        w_star = (n_l * y_l_bar + kappa * w)/(n_l + kappa)
        mu[l,] <- w_star + (chol(Sigma[,,l])/sqrt(kappa + n_l)) %*% rnorm(p)
      }
      else{
        ## sampling mean and variance
        y_l = y[z[s-1,] == l, ]
        n_l = sum(z[s-1,]==l)
        y_l_bar = colMeans(y_l)
        ## covariance matrix sampled
        V <- y_l - matrix(rep(y_l_bar,n_l), nrow = n_l, byrow = T)
        V <- crossprod(data.matrix(V))
        C_star <- C +  V + (n_l * kappa/(n_l + kappa)) * tcrossprod(y_l_bar - w)
        Sigma[,,l] <- rinvwishart(nu = n_l + r, S = C_star)
        Sigma_diag[l, ] <- diag(Sigma[,,l])
        Sigma_LT[l, ] <- Sigma[,,l][lower.tri(Sigma[,,l], diag = F)]
        ## mean sampled
        w_star = (n_l * y_l_bar + kappa * w)/(n_l + kappa)
        mu[l,] <- w_star + (chol(Sigma[,,l])/sqrt(kappa + n_l)) %*% rnorm(p)
      }
    }
    # sample cluster labels
    tau <- maketau(Pi = Pi, y = y, mu = mu, Sigma = Sigma)
    z[s,] <- apply(tau, 1, function(x) sample(L, size = 1, replace = T, prob = x))
    
    #for (i in 1:n){
      # for (l in 1:L) {
        # tau[l] <- log(Pi[l]) + mclust::dmvnorm(data=matrix(y[i,],nrow = 1), mean=matrix(mu[l,],nrow = 1), sigma=Sigma[,,l],log=T)
      #}
      # tau <- exp(tau - max(tau))
      # z[s,i] <- sample(1:L, size = 1, replace = T, prob = tau)
    #}
    # sub-sampling for theta samples
    z_mat <- model.matrix(~ -1 + cl, data = data.frame(cl = factor(z[s,], levels = 1:L)))
    total_pars <- cbind(mu, Sigma_diag, Sigma_LT)
    theta[[s]] <- z_mat %*% total_pars
    
    # print stops
    if (s %in% stops){
      print(s)
    }
  }
  return(list(theta = theta, z=z))
}