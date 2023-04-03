# evaluates the MERGELOSS function given Delta and psi
mergeloss <- function(psi, Delta_lt, w) {
  psi_matrix <- outer(psi, psi, function(x,y) as.integer(x==y))
  psi_ut <- psi_matrix[lower.tri(psi_matrix)]
  loss_members <- psi_ut * Delta_lt + w*(1-psi_ut) * (1-Delta_lt)
  loss <- sum(loss_members)
  return(loss)
}