clusterspirals <- function(N) {
  # N = c(300,200,100)
  # packages
  library(sn)
  library(ggplot2)
  library(dplyr)
  library(mcclust.ext)
  library(RColorBrewer)
  library(ggsci)
  library(kernlab)
  library(reshape2)
  library(KODAMA)
  # moons data
  d <- 2
  data <- KODAMA::spirals(n = c(N/3,N/3,N/3), sd = c(0.3,0.25,0.15))
  truth <- rep(1:3, each = N/3)
  y <- scale(data)
  K_0 <- 3
  # fitting a GMM
  source("rfuncts/mvnorm_gibbs.R")
  source("rfuncts/mnorm_D_apply.R")
  sourceCpp("rcppfuncts/mnorm_D_arma.cpp")
  S <- 15000 # iterations
  B <- 1000 # burnin
  L <- 30 # number of components
  sig_sq <- 0.5
  # fitting
  fit <- mvnorm_gibbs(S = S,
                      y = y,
                      L = L,
                      alpha = rep(1/2,L),
                      w = rep(0,d),
                      kappa = sig_sq/(2),
                      r = d+2,
                      C = sig_sq*diag(1,d),
                      stops = 50)
  # burn-out
  theta <- fit$theta[-(1:B)]
  M <- length(theta)
  n <- nrow(y)
  theta <- array(unlist(theta), dim = c(n, 2*d + choose(d,2), M))
  z <- fit$z[-(1:B),]
  # thin
  trip_ind <- seq(4,M,by=4)
  theta <- theta[,,trip_ind]
  z <- z[trip_ind,]
  # FOLD
  Delta_UT <- makeHellingerAvg(theta = theta, d = d, n = n)
  Delta <- matrix(0, nrow = n, ncol = n)
  Delta[upper.tri(Delta)] <- Delta_UT
  Delta <- Delta + t(Delta)
  Delta_H <- Delta
  # average-linkage
  Delta <- as.dist(Delta)
  cl <- hclust(Delta, method = "average")
  # list of possible clusterings
  max.k <- 10  # maximum number of clusters
  labs <- matrix(0, nrow = max.k, ncol = n)
  for (h in 1:max.k){
    labs[h,] = cutree(cl, k = h)
  }
  w_avg <- sum(Delta_UT)/(sum(1 - Delta_UT))
  losses <- apply(labs, 1, mergeloss, Delta_lt = Delta_H[lower.tri(Delta_H)], w=w_avg)
  # minimize
  min_loss <- which.min(losses)
  c.fold = labs[min_loss,]
  # credible ball and cGPSM
  cgsamps <- getcGsamps(theta=theta, w=w_avg, d=2, max.k=10)
  cb <- cGball(c.fold,cgsamps)
  # compare clusterings
  vi.fold <- mcclust::vi.dist(c.fold, truth)
  # evaluations
  inc <- vi.fold <= cb$dist.horiz
  k_fold <- length(table(c.fold))
  ari_fold <- adjustedRandIndex(c.fold, truth)
  logic_vec <- c(inc, k_fold, ari_fold)
  names(logic_vec) = c("Ball", "K_FOLD", "ARI")
  return(logic_vec)
}
