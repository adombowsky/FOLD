locnorm_D_apply <- function(theta, Sigma, d){
  mu <- theta
  D <- locnorm_D_arma(mu = mu,
                    Sigma = Sigma)
  D_UT <- t(D)[lower.tri(t(D))]
  return(D_UT)
}