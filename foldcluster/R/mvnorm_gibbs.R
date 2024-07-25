#' Gibbs sampler for a multivariate Bayesian GMM
#'
#' Samples from the posterior distribution of a GMM with unknown means and covariances.
#' @param S Total number of MCMC iterations.
#' @param y Matrix with \code{n} observations as the rows and \code{p} features as the columns.
#' @param L Number of mixture components. Default is \code{L=10}.
#' @param alpha Vector of Dirichlet concentration parameters. Default is \code{(1/L, ..., 1/L)}.
#' @param w Normal-Inverse-Wishart location parameter, \code{p}-dimensional.
#' @param kappa Normal-Inverse-Wishart scale parameter.
#' @param r Inverse-Wishart degrees of freedom, greater than \code{p}-1.
#' @param C Inverse-Wishart scale matrix, \code{p}x\code{p}.
#' @param stops Number of iterations to print progress point.
#' @importFrom Rcpp sourceCpp
#' @importFrom LaplacesDemon rdirichlet rinvwishart
#' @importFrom stats rnorm
#' @return A list MCMC samples from the posteriors of: \code{theta} (localized atoms), and \code{z} (component labels).
#' @details This function fits a Bayesian GMM in order to cluster a sample of observations. The prior distributions are
#' conjugate; the means and covariance are drawn from a normal-Inverse-Wishart distribution with parameters \code{w},
#' \code{kappa}, \code{r}, and \code{C}, while the weights follow a Dirichlet distribution with parameter \code{alpha}.
#' @export
#' @useDynLib foldcluster
mvnorm_gibbs <- function(S, y, L=10, alpha=rep(1/L,L), w, kappa, r, C, stops = 1){

  # prelims
  n <- nrow(y)
  p <- ncol(y)
  # initialization
  Pi <- alpha/sum(alpha) # weights
  mu <- matrix(0, nrow = L, ncol = p) # component means
  Sigma <- array(diag(1,p), dim = c(p,p,L))
  Sigma_diag <- matrix(1, nrow = L, ncol = p) # component variances
  Sigma_LT <- matrix(1, nrow = L, ncol = choose(p,2)) # lower triangles
  z <- matrix(sample(1:L, size = n*S, replace = T), nrow = S, ncol = n ) # cluster labels
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
        C_star <- C +  (n_l * kappa/(n_l + kappa)) * tcrossprod(y_l_bar - w)
        Sigma[,,l] <- LaplacesDemon::rinvwishart(nu = n_l + r, S = C_star)
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
        Sigma[,,l] <- LaplacesDemon::rinvwishart(nu = n_l + r, S = C_star)
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

    # save output
    total_pars <- cbind(mu, Sigma_diag, Sigma_LT)
    theta[[s]] <- total_pars[z[s,],]

    # print stops
    if (s %in% stops){
      print(s)
    }
  }
  return(list(theta = theta, z=z))
}

#' Gibbs sampler for a univariate Bayesian GMM
#'
#' Samples from the posterior distribution of a GMM with unknown means and variances.
#' @param S Total number of MCMC iterations.
#' @param y Vector with \code{n} observations to cluster.
#' @param L Number of mixture components. Default is \code{L=10}.
#' @param alpha Vector of Dirichlet concentration parameters. Default is \code{(1/L, ..., 1/L)}.
#' @param w Normal-Inverse-Gamma location parameter, scalar.
#' @param kappa Normal-Inverse-Gamma scale parameter, positive.
#' @param r Inverse-Gamma degrees of freedom, positive.
#' @param C Inverse-Gamma scale, positive.
#' @param stops Number of iterations to print progress point.
#' @importFrom Rcpp sourceCpp
#' @importFrom LaplacesDemon rdirichlet rinvwishart
#' @importFrom stats rnorm
#' @return A list MCMC samples from the posteriors of: \code{theta} (localized atoms), and \code{z} (component labels).
#' @details This function fits a univariate Bayesian GMM in order to cluståer a sample of observations. The prior distributions are
#' conjugate; the means and variances are drawn from a normal-Inverse-Gamma distribution with parameters \code{w},
#' \code{kappa}, \code{r}, and \code{C}. Note that we have parameterized the variance prior as
#' Inv-Gamma(\code{r}/2,\code{C}/2), in order to be consistent with the normal-Inverse-Wishart for higher dimensions.
#' The weights follow a Dirichlet distribution with parameter \code{alpha}.
#' @export
#' @useDynLib foldcluster
unorm_gibbs <- function(S, y, L=10, alpha=rep(1/L,L), w, kappa, r, C, stops = 1){

  # prelims
  n <- length(y)
  p <- 1
  # initialization
  Pi <- alpha/sum(alpha) # weights
  mu <- rep(0,L) # component means
  Sigma <- rep(1,L) # component variances
  z <- matrix(sample(1:L, size = n*S, replace = T), nrow = S, ncol = n ) # cluster labels
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
        Sigma[l] <- LaplacesDemon::rinvwishart(nu = r, S = C)
        mu[l] <- w + (chol(Sigma[l])/sqrt(kappa)) %*% rnorm(p)
      }
      else if (cluster.sizes[l] == 1) {
        ## sampling mean and covariance
        y_l = y[z[s-1,] == l]
        n_l = sum(z[s-1,]==l)
        y_l_bar = y_l
        ## covariance matrix sampled
        C_star <- C +  (n_l * kappa/(n_l + kappa)) * tcrossprod(y_l_bar - w)
        Sigma[l] <- LaplacesDemon::rinvwishart(nu = n_l + r, S = C_star)
        ## mean sampled
        w_star = (n_l * y_l_bar + kappa * w)/(n_l + kappa)
        mu[l] <- w_star + (chol(Sigma[l])/sqrt(kappa + n_l)) %*% rnorm(p)
      }
      else{
        ## sampling mean and variance
        y_l = y[z[s-1,] == l]
        n_l = sum(z[s-1,]==l)
        y_l_bar = mean(y_l)
        ## covariance matrix sampled
        V <- y_l - matrix(rep(y_l_bar,n_l), nrow = n_l, byrow = T)
        V <- crossprod(data.matrix(V))
        C_star <- C +  V + (n_l * kappa/(n_l + kappa)) * tcrossprod(y_l_bar - w)
        Sigma[l] <- LaplacesDemon::rinvwishart(nu = n_l + r, S = C_star)
        ## mean sampled
        w_star = (n_l * y_l_bar + kappa * w)/(n_l + kappa)
        mu[l] <- w_star + (chol(Sigma[l])/sqrt(kappa + n_l)) %*% rnorm(p)
      }
    }
    # sample cluster labels
    tau <- makeutau(Pi = Pi, y = y, mu = mu, Sigma = Sigma)
    z[s,] <- apply(tau, 1, function(x) sample(L, size = 1, replace = T, prob = x))

    # save output
    total_pars <- cbind(mu, Sigma)
    theta[[s]] <- total_pars[z[s,],]

    # print stops
    if (s %in% stops){
      print(s)
    }
  }
  return(list(theta = theta, z=z))
}

#' Gibbs sampler for a location Bayesian GMM with fixed covariance
#'
#' Samples from the posterior distribution of a location GMM where the covariance is fixed and constant.
#' @param S Total number of MCMC iterations.
#' @param y Matrix with \code{n} observations as the rows and \code{p} features as the columns.
#' @param L Number of mixture components. Default is \code{L=10}.
#' @param alpha Vector of Dirichlet concentration parameters. Default is \code{(1/L, ..., 1/L)}.
#' @param mu_0 Prior mean for component centers.
#' @param Sigma_0 Prior covariance for component centers.
#' @param Sigma Covariance of the data, fixed and constant across components.
#' @param stops Number of iterations to print progress point.
#' @importFrom Rcpp sourceCpp
#' @importFrom LaplacesDemon rdirichlet
#' @importFrom stats rnorm
#' @return A list MCMC samples from the posteriors of: \code{theta} (localized centers), and \code{z} (component labels).
#' @details This function fits a special case of the Bayesian GMM in which the component-specific covariance \code{Sigma} is
#' fixed and constant across all components. The cluster means, or centers, are drawn from a normal distribution with parameters
#' \code{mu_0} and \code{Sigma_0}, while the weights are given a Dirichlet prior with parameter \code{alpha}.
#' @export
#' @useDynLib foldcluster
lmvnorm_gibbs <- function(S, y, L, alpha=rep(1/L,L), mu_0, Sigma_0, Sigma, stops = 1){
  # prelims
  n <- nrow(y)
  p <- ncol(y)
  # initialization
  Pi <- alpha/sum(alpha) # weights
  mu <- matrix(0, nrow = L, ncol = p) # component means
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
    tau <- lmaketau(Pi = Pi, y = y, mu = mu, Sigma = Sigma)
    z[s,] <- apply(tau, 1, function(x) sample(L, size = 1, replace = T, prob = x))
    # sub-sampling for theta samples
    theta[[s]] <- cbind(mu[z[s,],],
                        do.call(cbind,as.list(diag(Sigma)))[rep(1,n),],
                        do.call(cbind,as.list(rep(0,choose(p,2))))[rep(1,n),])

    # print stops
    if (s %in% stops){
      print(s)
    }
  }
  return(list(theta = theta, z=z))
}
