# computes the squared Hellinger distance between P~N(mu_1, Sigma_1), Q~N(mu_2, Sigma_2)
mnorm.hellinger <- function(mu_1, Sigma_1, mu_2, Sigma_2){
  C <- det(Sigma_1)^(1/4) * det(Sigma_2)^(1/4) / (det( (Sigma_1 + Sigma_2)/2   )^(1/2))
  M <- -(1/8) * crossprod(mu_1 - mu_2, solve((Sigma_1 + Sigma_2)/2) %*% (mu_1 - mu_2))
  H <- 1 - C * exp(M)
  return(H)
}