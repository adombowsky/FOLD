nongauss.cluster <- function(n, S, B, d, L){
  Rcpp::sourceCpp("rcppfuncts/mnorm_D_arma.cpp")
  ari <- c()
  K <- c()
  samp <- nongauss(n = n,
                          Pi = c(0.2, 0.3, 0.15, 0.35),
                          xi = matrix(c(2.5,3.5,
                                        0,-3.5), byrow = T, nrow = 2),
                          Omega = array(c(1,0,0,1,5,0,0,2), c(2,2,2)),
                          alpha = matrix(c(-10,15,
                                           4,-17), byrow = T, nrow = 2),
                          Mu_1 = c(-4,2.5),
                          Sigma_1 = matrix(c(1/2, 1/2, 1/2, 5/2), byrow = T, nrow = 2),
                          # Mu = matrix(rnorm(4, mean = 1, sd = 1), nrow = 2),
                          Mu = matrix(c(2.325, 1.085, 4.381, 2.009), nrow = 2),
                          Sigma = array(c(0.2, 0, 0, 0.8, 0.7, 0, 0, 0.6), c(2,2,2)),
                          p = 1/3)
  s <- samp$s # true labels
  # relabeling s
  s[s==4] = rep(1, sum(s==4))
  y <- samp$data # data
  y <- scale(y) # scaled
  # fitting
  fit <- mvnorm_gibbs(S = S,
                      y = y,
                      L = L,
                      alpha = rep(1/2,L), # rep((d/2)+2,L)
                      w = rep(0,d),
                      kappa = 1,
                      r = d+2,
                      C = diag(1,d),
                      stops = S)
  # burn-out
  theta <- fit$theta[-(1:B)]
  M <- length(theta)
  n <- nrow(y)
  theta <- array(unlist(theta), dim = c(n, 2*d + choose(d,2), M))
  z <- fit$z[-(1:B),]
  # thin
  trip_ind <- seq(3,M,by=3)
  theta <- theta[,,trip_ind]
  z <- z[trip_ind,]
  # apply hellinger distance function
  Delta_UT <- makeHellingerAvg(theta = theta, d = d, n = n)
  Delta <- matrix(0, nrow = n, ncol = n)
  Delta[upper.tri(Delta)] <- Delta_UT
  Delta <- Delta + t(Delta)
  Delta_H <- Delta
  # complete-linkage
  Delta <- as.dist(Delta)
  cl <- hclust(Delta, method = "complete")
  # list of possible clusterings
  max.k <- 20  # maximum number of clusters
  labs <- matrix(0, nrow = max.k, ncol = n)
  for (h in 1:max.k){
    labs[h,] = cutree(cl, k = h)
  }
  w_avg <- sum(Delta_UT)/(sum(1 - Delta_UT))
  gamma_med <- median(Delta_UT)
  w_med = gamma_med/(1-gamma_med)
  losses <- apply(labs, 1, mergeloss, Delta_lt = Delta_H[lower.tri(Delta_H)], w=w_avg)
  # minimize
  min_loss <- which.min(losses)
  c.fold = labs[min_loss,]
  K[1] <- length(unique(c.fold))
  ari[1] <- round(adjustedRandIndex(c.fold,s),3)
  # W & G estimate
  c.psm <- mcclust::comp.psm(z)
  c.mv <- mcclust.ext::minVI(psm = c.psm,
                             cls.draw = z,
                             max.k = 20)
  c.VI <- c.mv$cl
  K[2] <- length(unique(c.VI))
  ari[2] <- round(adjustedRandIndex(c.VI,s),3)
  # W & G estimate -- Binder's Loss
  c.b <- mcclust::minbinder(psm = c.psm,
                            cls.draw = z,
                            max.k = 20)
  c.Binder <- c.b$cl 
  K[3] <- length(unique(c.Binder))
  ari[3] <- round(adjustedRandIndex(c.Binder,s),3)
  # Mclust
  mcl <- mclust::Mclust(y, G=1:max.k)
  c.mcl <- mcl$classification
  K[4] <- length(unique(c.mcl))
  ari[4] <- round(adjustedRandIndex(c.mcl,s),3)
  return(list(K=K,ari=ari))
  return(list(K=K,ari=ari))
}