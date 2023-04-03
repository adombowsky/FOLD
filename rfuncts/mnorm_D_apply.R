mnorm_D_apply <- function(theta, d){
  mu <- theta[ , 1:d]
  Sigma_diag <- theta[ , ((d+1):(2*d))]
  Sigma_LT <- theta[, -(1:(2*d))]
  if (d == 2) {
    Sigma_LT <- matrix(Sigma_LT, ncol = 1)
  }
  D <- mnorm_D_arma(mu = mu,
                    Sig_diag = Sigma_diag,
                    Sig_LT = Sigma_LT)
  D_UT <- t(D)[lower.tri(t(D))]
  return(D_UT)
}