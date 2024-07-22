#' @importFrom Rcpp sourceCpp
#' @importFrom stats hclust cutree
makecG <- function(theta, omega, p, max.k){
  if (p==1) {
    mu <- theta[,1]
    sigma <- sqrt(theta[,2])
    D <- unorm_D_arma(mu=mu, sigma=sigma)
  } else {
    mu <- theta[,1:p]
    Sig_diag <- theta[, ((p+1):(2*p))]
    Sig_LT <- theta[,-(1:(2*p))]
    if (p==2) {
      Sig_LT <- as.matrix(Sig_LT)
    }
    D <- mnorm_D_arma(mu=mu, Sig_diag = Sig_diag, Sig_LT = Sig_LT)
  }
  D <- D + t(D)
  D.dist <- as.dist(D)
  cl <- hclust(D.dist, method = "average")
  labs <- matrix(0, nrow = max.k, ncol = n)
  for (h in 1:max.k){
    labs[h,] = cutree(cl, k = h)
  }
  c.G <- min_risk(c=labs, Delta = D, omega=omega)
  return(c.G)
}

#' Posterior samples of \code{c_theta}.
#'
#' @param theta A list of MCMC samples from the posterior of the localized atoms.
#' @param omega The FOLD separation parameter, positive.
#' @param p Number of features.
#' @param max.k Maximum number of clusters.
#' @return A matrix with rows corresponding to samples from the posterior of \code{c_theta}.
#' @details This function extracts posterior samples from \code{c_theta}, a minimizer of the FOLD
#'  loss function. The resulting matrix of samples can then be used to compute the credible ball.
#' @export
#' @useDynLib fold
get_ctheta_samps <- function(theta, omega, p, max.k) {
  theta <- array(unlist(theta), dim = c(n, 2*p + choose(p,2), length(theta)))
  cgsamps <- t(apply(theta, 3, function(x) makecG(x, omega=omega, p=p, max.k=max.k)))
  return(cgsamps)
}

#' Credible ball for a FOLD clustering point estimate.
#'
#' @param c_fold A FOLD clustering point estimate.
#' @param ctheta_samps A matrix of posterior samples for \code{c_theta}.
#' @return A list containing the bounds of the 95% credible ball around \code{c_fold}.
#' @details Using a point estimate and the posterior samples for \code{c_theta}, this function
#'  constructs a 95% credible ball. The credible ball can be summarized by the horizontal bounds
#'  (\code{c.horiz}), the upper vertical bounds (\code{c.uppervert}), or the lower vertical bounds
#'  (\code{c.lowervert}).
#' @export
#' @useDynLib fold
ctheta_ball <- function(c_fold, ctheta_samps){
  cb <- credibleball(c_fold,ctheta_samps)
  return(cb)
}
