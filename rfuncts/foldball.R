makecG <- function(theta, w, d, max.k){
  require(mcclust.ext)
  require(Rcpp)
  require(RcppArmadillo)
  source("rfuncts/mergeloss.R")
  mu <- theta[,1:d]
  Sig_diag <- theta[, ((d+1):(2*d))]
  Sig_LT <- theta[,-(1:(2*d))]
  if (d==2) {
    Sig_LT <- as.matrix(Sig_LT)
  }
  D <- mnorm_D_arma(mu=mu, Sig_diag = Sig_diag, Sig_LT = Sig_LT)
  D <- D + t(D)
  D.dist <- as.dist(D)
  cl <- hclust(D.dist, method = "average")
  labs <- matrix(0, nrow = max.k, ncol = n)
  for (h in 1:max.k){
    labs[h,] = cutree(cl, k = h)
  }
  losses <- apply(labs, 1, mergeloss, Delta_lt = Delta_H[lower.tri(Delta_H)], w=w)
  # minimize
  min_loss <- which.min(losses)
  c.G <- labs[min_loss,]
  return(c.G)
}

getcGsamps <- function(theta, w, d, max.k) {
  cgsamps <- t(apply(theta, 3, function(x) makecG(x, w=w, d=d, max.k=max.k)))
  return(cgsamps)
}

cGball <- function(c.fold, cgsamps){
  require(mcclust)
  require(mcclust.ext)
  # credible ball
  cb <- credibleball(c.fold,cgsamps)
  return(cb)
}

cGpsm <- function(cgsamps){
  require(mcclust)
  require(mcclust.ext)
  # credible ball
  cb <- mcclust::comp.psm(cgsamps)
  return(cb)
}
