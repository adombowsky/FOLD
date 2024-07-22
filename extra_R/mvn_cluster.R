mvn.cluster <- function(n, S, B, d, L){
  Rcpp::sourceCpp("extra_src/oracle.cpp")
  ari <- c()
  K <- c()
  mu <- matrix(c(6.5,5,
                 0,0,
                 -5,-5), byrow = T, nrow = 3)
  Sigma <- array(c(1,0,0,1,5,0,0,2,3,0,0,1), c(2,2,3))
  Pi <- c(0.45, 0.25, 0.3)
  samp <- multivar.normal(n =  n,
                          Pi = Pi,
                          xi = mu,
                          Omega = Sigma)
  s <- samp$s # true labels
  y_us <- samp$data # data
  y <- scale(samp$data) # scaled
  # oracle clustering parameters
  h <- makeh(mu, Sigma)
  Deltastar <- makeOracle(X=y_us, mu = mu, Sigma = Sigma, pi = Pi, h = h)
  Deltastar <- Deltastar + t(Deltastar)
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
  n <- nrow(y)
  z <- fit$z[-(1:B),]
  # thin
  trip_ind <- seq(3,nrow(z),by=3)
  theta <- theta[trip_ind]
  z <- z[trip_ind,]
  # apply hellinger distance function
  Delta <- comp_delta(theta,d,n)
  # complete-linkage
  cl <- hclust(as.dist(Delta), method = "complete")
  clstar <- hclust(as.dist(Deltastar), method = "complete")
  # list of possible clusterings
  max.k <- 20  # maximum number of clusters
  labs <- labs_star <- matrix(0, nrow = max.k, ncol = n)
  for (h in 1:max.k){
    labs[h,] = cutree(cl, k = h)
    labs_star[h,] = cutree(clstar, k = h)
  }
  Delta_UT <- Delta[upper.tri(Delta)]
  w_avg <- sum(Delta_UT)/(sum(1 - Delta_UT))
  c.fold <- min_risk(labs, Delta, w_avg)
  c.star <- min_risk(labs_star, Deltastar, w_avg)
  K[1] <- length(unique(c.fold))
  ari[1] <- round(adjustedRandIndex(c.fold,s),3)
  K[2] <- length(unique(c.star))
  ari[2] <- round(adjustedRandIndex(c.star,s),3)
  # W & G estimate
  c.psm <- mcclust::comp.psm(z)
  c.mv <- mcclust.ext::minVI(psm = c.psm,
                             cls.draw = z,
                             max.k = 20)
  c.VI <- c.mv$cl
  K[3] <- length(unique(c.VI))
  ari[3] <- round(adjustedRandIndex(c.VI,s),3)
  # W & G estimate -- Binder's Loss
  c.b <- mcclust::minbinder(psm = c.psm,
                            cls.draw = z,
                            max.k = 20)
  c.Binder <- c.b$cl
  K[4] <- length(unique(c.Binder))
  ari[4] <- round(adjustedRandIndex(c.Binder,s),3)
  # Mclust
  mcl <- mclust::Mclust(y, G=1:max.k)
  c.mcl <- mcl$classification
  K[5] <- length(unique(c.mcl))
  ari[5] <- round(adjustedRandIndex(c.mcl,s),3)
  return(list(K=K,ari=ari))
}
