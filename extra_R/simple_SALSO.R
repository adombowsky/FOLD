r_simple_SALSO <- function(Delta, omega, p_SA=0.5, max.k=10, max.zealous=10) {
  # RCPP
  require(Rcpp)
  require(RcppArmadillo)
  sourceCpp("extra_src/salso.cpp")
  n <- nrow(Delta)
  # initialization
  u <- runif(1)
  if (u < p_SA) {
    ind0 = sample(1:n)
    c0 = rep(0,n)
    K0 = 1
    c0[1] = 1
    risk0 = c()
    minrisk0 = 0
    for (i in 2:n) {
      Delta_i = Delta[ind0[1:i], ind0[1:i]]
      for (k in 1:(K0+1)) {
        c0[i] = k
        risk0[k] = risk_cpp(c0[1:i], Delta_i, omega)
      }
      minrisk0 = which.min(risk0)
      c0[i] = minrisk0
      if (minrisk0 == (K0+1)) {
        K0 = K0+1
      }
    }
    c0 = c0[order(ind0)]
  } else {
    c0 = sample(max.k, n, replace = TRUE)
    K0 = max.k
  }
  # sweetening
  c1 = c0
  ind1 = sample(1:n)
  K1 = K0
  risk1 = c()
  minrisk1 = 0
  Delta1 = Delta[ind1,ind1]
  for (i in 1:n) {
    for (k in 1:(K1+1)) {
        c1[i] = k
        risk1[k] = risk_cpp(c1, Delta1, omega)
    }
      minrisk1 = which.min(risk1)
      c1[i] = minrisk1
      if (minrisk1 == (K1+1)) {
        K1 = K1+1
    }
  }
  c1 = c1[order(ind1)] # re-ordering
  c1 = as.numeric(as.factor(c1)) # re-normalizing
  K1 = max(c1)
  # if (rand_index(c1,c0)==1) { break }
  # zealous updates
  K2 = K1
  cl_ind = sample(1:K2)
  c2 = c1
  n_out = 0
  risk2 = c()
  minrisk2 = 0
  for (m in 1:min(max.zealous, K1)) {
    c2_old = c2
    c2_out = which(c2==cl_ind[m])
    # re-allocation
    ind_out = sample(c2_out)
    for (i in ind_out) {
      for (k in 1:(K2+1)) {
        c2[i]  = k
        risk2[k] = risk_cpp(c2, Delta, omega)
      }
      minrisk2 = which.min(risk2)
      c2[i] = minrisk2
      if (minrisk2 == (K2+1)) {
        K2 = K2+1
      }
    }
    if (risk_cpp(c2,Delta,omega) >= risk_cpp(c2_old, Delta, omega)) {
      c2 = c2_old
    }
  }
  c2 = as.numeric(as.factor(c2))
  return(c2)
}
