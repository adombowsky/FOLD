convergence_line <- function(n, seed,S){
  # packages and functions
  require(fold)
  require(Rcpp)
  require(RcppArmadillo)
  require(ggplot2)
  require(cowplot)
  source("extra_R/simple_SALSO.R")
  source("extra_R/salso_parallel.R")
  source("extra_R/locationnorm_fixed_sigma_gibbs.R")
  sourceCpp("extra_src/locnorm_rsig_D_arma.cpp")
  sourceCpp("extra_src/oracle_convergence.cpp")
  sourceCpp("extra_src/salso.cpp")
  # step 1: simulate data
  set.seed(1)
  d <- 2
  K0 <- 4
  mu0 <- matrix(c(1, 1,
                  1.75, 1.75,
                  -1.75, -1.75,
                  -1,-1),
                ncol = d,
                byrow=TRUE)
  x <- matrix(0, nrow = n, ncol = d)
  Pi0 <- rep(1/K0,K0)
  s0 <- sample(K0, size = n, replace = TRUE)
  tau <- 0.5
  for (i in 1:n) {
    x[i,] <- mu0[s0[i],] + tau*rnorm(d)
  }
  B <- S/2 # burnin
  L <- 4
  # fitting
  fit <- loc_fsig_mvnorm_gibbs_cpp(S = S,
                                   y = x,
                                   L = L,
                                   alpha = rep(1/L,L),
                                   mu_0 = rep(0,d),
                                   Sigma_0 = diag(1,d),
                                   sigma_sq = tau^2,
                                   stops = S)
  theta <- fit$theta[-(1:B)]
  sigma_sq <- tau^2
  M <- length(theta)
  theta <- array(unlist(theta), dim = c(n, d, M))
  z <- fit$z[-(1:B),]
  # thin
  trip_ind <- seq(2,M,by=2)
  theta <- theta[,,trip_ind]
  z <- z[trip_ind,]
  # apply hellinger distance function
  Delta_UT <- locrsigHellingerAvg(theta = theta, sigma_sq = sigma_sq*rep(1,dim(theta)[3]), d = d, n = n)
  Delta <- matrix(0, nrow = n, ncol = n)
  Delta[upper.tri(Delta)] <- Delta_UT
  Delta <- Delta + t(Delta)
  Delta_H <- Delta
  # compute oracle values
  Sigma0 <- array(tau^2*diag(d), dim = c(d,d,K0))
  h <- makeh(mu0, Sigma0)
  Deltastar <- makeOracle(X=x, mu = mu0, Sigma = Sigma0, pi = Pi0, h = h)
  Deltastar_UT <- Deltastar
  Deltastar <- Deltastar + t(Deltastar)
  Deltastar_H <- Deltastar
  # average-linkage
  Delta <- as.dist(Delta)
  Deltastar <- as.dist(Deltastar)
  cl <- hclust(Delta, method = "average")
  clstar <- hclust(Deltastar, method = "average")
  # compute clusterings
  max.k <- 10  # maximum number of clusters
  labs <- labs_star <- matrix(0, nrow = max.k, ncol = n)
  for (h in 1:max.k){
    labs[h,] = cutree(cl, k = h)
    labs_star[h,] = cutree(clstar, k = h)
  }
  w_avg <- sum(Delta_UT)/(sum(1 - Delta_UT))
  w_avg_star <- sum(Deltastar_UT[upper.tri(Deltastar_UT)])/(sum(1-Deltastar_UT[upper.tri(Deltastar_UT)]))
  losses <- apply(labs, 1, mergeloss, Delta_lt = Delta_H[lower.tri(Delta_H)], w=w_avg)/choose(n,2)
  losses_star <- apply(labs, 1, mergeloss, Delta_lt = Deltastar_H[lower.tri(Deltastar_H)], w = w_avg_star)/choose(n,2)
  # minimize
  line.df <- data.frame(K=rep(1:max.k,2),
                        loss = c(losses, losses_star),
                        Risk = c(rep("FOLD", max.k), rep("Oracle", max.k)))
  line.plot <- ggplot(line.df, aes(x = K, y = loss, col = Risk)) + geom_line() + geom_point() +
    ylab("Risk") +
    xlab("Number of Clusters") +
    scale_x_discrete(limits = as.character(1:max.k)) +
    labs(title = paste0("n=", toString(n))) +
    theme_bw() +
    theme(text = element_text(size=15))
  # computing w/ salso
  c.salso = salso_parallel(M=20,
                           Delta=Delta_H,
                           omega = w_avg)
  c.salso.oracle = salso_parallel(M=20,
                                  Delta=Deltastar_H,
                                  omega = w_avg_star)
  # compare to Binder
  Beta = mcclust::comp.psm(z)
  c.B = salso_parallel(M=20,
                       Delta = 1-Beta,
                       omega = 1)
  # compare to Binder oracle
  Betastar <- makeBinderOracle(X=x, mu = mu0, Sigma = Sigma0, pi = Pi0)
  Betastar_UT <- Betastar
  Betastar <- Betastar + t(Betastar)
  c.B.oracle = salso_parallel(M=20,
                              Delta=Betastar,
                              omega = 1)
  return(list(plot = line.plot,
         data = x,
         c.fold = c.salso,
         c.oracle = c.salso.oracle,
         c.B = c.B,
         c.B.oracle = c.B.oracle))
}
