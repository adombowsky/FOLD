salso_parallel <- function(M, Delta, omega, p_SA=0.5, max.k=10, max.zealous = 10) {
  ## packages
  library(foreach)
  library(doParallel)
  library(Rcpp)
  library(RcppArmadillo)
  source("extra_R/simple_SALSO.R")
  ## configuring
  n.cores <- 10
  ### create the cluster
  my.cluster <- parallel::makeCluster(
    n.cores,
    type = "PSOCK"
  )
  ### register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  x <- foreach(r = 1:M,
            .export = c("r_simple_SALSO"),
            .noexport = c("risk_cpp",
                          "rand_index")) %dopar% {
                            r_simple_SALSO(Delta=Delta, omega=omega, p_SA = p_SA, max.k = max.k, max.zealous = max.zealous)
                          }
  partitions = matrix(unlist(x), nrow = M, byrow = TRUE)
  risks <- risk_matrix_cpp(partitions, Delta, omega)
  c.star = partitions[which.min(risks),]
  return(c.star)
}
