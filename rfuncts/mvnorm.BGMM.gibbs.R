gmm.gibbs <- function(R, B, y, alpha, kappa, w, nu, M, K, stops = 1) {
  # Gibbs sampler from McLachlan and Peel
  require(mvtnorm)
  require(LaplacesDemon)
  
  # K = number of clusters
  # alpha = Dirichlet parameter
  # kappa = scaling parameter for GMM prior
  # w = matrix of mean parameters for mu's
  # nu = first parameter in Wishart distribution
  # M = second parameter in Wishart distribution
  
  # prelims
  n <- nrow(y)
  p <- ncol(y)
  # storage
  pi <- matrix(NA, nrow = R, ncol = K)
  mu <- list() # list of matricies for each cluster
  for (i in 1:R) {
    mu[[i]] <- matrix(0, nrow = K, ncol = p)
  }
  Sigma <- list() # list of list, each inner list containing Sigma_k samples
  for (r in 1:R) {
    Sigma[[r]] <- list() 
  }
  c <- matrix(NA, nrow = R, ncol = n) # cluster labels
  cluster.means <- matrix(NA, nrow = K, ncol = p)
  # initialization
  Z <- t(rmultinom(n, 1, rep(1/K, K)))
  c[1, ] <- apply(Z, 1, function(x) which.max(x == 1))
  M_inv <- solve(M)
  tau <- c()
  for (k in 1:K){
    Sigma[[1]][[k]] <- diag(p)
  }
  # configuring stops
  stops <- (1:(R/stops)) * stops
  # sampling
  print("Sampling")
  for (r in 2:R) {
    # compute cluster sizes
    cluster.sizes <- as.vector(table(factor(c[r-1, ], levels = 1:L)))
    # sample pi
    alpha_star <- alpha + cluster.sizes
    pi[r, ] <- LaplacesDemon::rdirichlet(1, alpha_star)
    # sample normal parameters
    for (k in 1:K){
      # mu
      cluster <- as.matrix(y[c[r-1, ] == k, ])
      cluster.means[k, ] <- ifelse(cluster.sizes[k] > 0, colMeans(cluster), rep(0,p))
      w_star = (cluster.sizes[k] * cluster.means[k, ] + kappa * w)/(cluster.sizes[k] + kappa)
      mu[[r]][k, ] <- rmvnorm(1, w_star, (1/(cluster.sizes[k] + kappa)) * Sigma[[r-1]][[k]])
      # Sigma
      V <- ifelse(cluster.sizes[k]>0,
                  tcrossprod(cluster[1, ] - cluster.means[k, ], cluster[1, ] - cluster.means[k, ]),
                  matrix(0, nrow = p, ncol = p))
      if (cluster.sizes[k] > 1) {
        for (j in 2:cluster.sizes[k]){
          V <- V + tcrossprod(cluster[j, ] - cluster.means[k, ], cluster[j, ] - cluster.means[k, ])
        }
      }
      M_star = M_inv +  V + ( (cluster.sizes[k] * nu)/(cluster.sizes[k] + nu)) * 
        tcrossprod(cluster.means[k, ] - w, cluster.means[k, ] - w)
      M_star = round(M_star, 6)
      Sigma[[r]][[k]] <- rinvwishart(cluster.sizes[k] + nu, M_star)
    }
    # sample multinomial vectors and cluster labels
    for (i in 1:n) {
      for (k in 1:K) {
        tau[k] <- pi[r, k] * dmvnorm(x = y[i, ], mean = mu[[r]][k, ], sigma = Sigma[[r]][[k]])
      }
      tau <- tau/sum(tau)
      Z[i, ] <- rmultinom(1,1,tau)
    }
    c[r, ] <- apply(Z, 1, function(x) which(x == 1))
    
    # print stops
    if (r %in% stops){
      print(r)
    }
  }
  
  # discard Burn-in
  pi <- pi[(B+1):R, ]
  mu <- mu[(B+1):R]
  Sigma <- Sigma[(B+1):R]
  c <- c[(B+1):R, ]
  return(list(pi = pi, mu = mu, Sigma = Sigma, c = c))
}