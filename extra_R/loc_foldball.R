makecG_loc <- function(theta, Sigma, w, d, max.k){
  require(foldcluster)
  require(Rcpp)
  require(RcppArmadillo)
  Rcpp::sourceCpp("extra_src/")
  D <- locnorm_D_arma(mu=theta, Sigma=Sigma)
  D <- D + t(D)
  cl <- hclust(as.dist(D), method = "average")
  labs <- matrix(0, nrow = max.k, ncol = n)
  for (h in 1:max.k){
    labs[h,] = cutree(cl, k = h)
  }
  c.G <- min_risk(labs, D, w)
  return(c.G)
}

getcGsamps_loc <- function(theta, Sigma, w, d, max.k) {
  cgsamps <- t(apply(theta, 3, function(x) makecG_loc(x, Sigma=Sigma, w=w, d=d, max.k=max.k)))
  return(cgsamps)
}
