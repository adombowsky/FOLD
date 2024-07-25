clusterspirals <- function(N) {
  # packages
  library(foldcluster)
  library(mclust)
  library(mcclust.ext)
  library(KODAMA)
  # moons data
  d <- 2
  data <- KODAMA::spirals(n = c(N/3,N/3,N/3), sd = c(0.3,0.25,0.15))
  truth <- rep(1:3, each = N/3)
  y <- scale(data)
  n <- nrow(y)
  K_0 <- 3
  # fitting a GMM
  S <- 15000 # iterations
  B <- 1000 # burnin
  L <- 30
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
  z <- fit$z[-(1:B),]
  # thin
  trip_ind <- seq(4,nrow(z),by=4)
  theta <- theta[trip_ind]
  z <- z[trip_ind,]
  # FOLD
  Delta <- comp_delta(theta,d,n)
  cl <- hclust(as.dist(Delta), method = "average")
  # list of possible clusterings
  max.k <- 10  # maximum number of clusters
  labs <- matrix(0, nrow = max.k, ncol = n)
  for (h in 1:max.k){
    labs[h,] = cutree(cl, k = h)
  }
  Delta_UT <- Delta[upper.tri(Delta)]
  w_avg <- sum(Delta_UT)/(sum(1 - Delta_UT))
  c.fold <- min_risk(labs,Delta,w_avg)
  # compute credible ball
  ctheta_samps <- get_ctheta_samps(theta=theta, omega=w_avg, p=d, max.k=max.k)
  cb <- ctheta_ball(c.fold,ctheta_samps)
  # compare clusterings
  vi.fold <- mcclust::vi.dist(c.fold, truth)
  k_horiz <- length(table(cb$c.horiz[1,]))
  k_uv <- length(table(cb$c.uppervert[1,]))
  k_lv <- length(table(cb$c.lowervert[1,]))
  # evaluations
  inc <- vi.fold <= cb$dist.horiz
  k_int <- (K_0 <= k_lv) & (K_0 >= k_uv)
  k_fold <- length(table(c.fold))
  ari_fold <- adjustedRandIndex(c.fold, truth)
  logic_vec <- c(inc, k_int, k_fold, ari_fold)
  names(logic_vec) = c("Ball", "K_interval", "K_FOLD", "ARI")
  return(logic_vec)
}
