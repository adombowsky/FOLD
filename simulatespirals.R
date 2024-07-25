# packages
library(foldcluster)
library(grid)
library(gridExtra)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doParallel)
# r functions
source("extra_R/clusterspirals.R")
# parallel computing
n.cores <- 10
## create the cluster
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
## register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
# parallel loop
N <- 300
R <- 100
x <- foreach(r = 1:R,
             .packages = c("mcclust", "Rcpp", "cluster", "KODAMA", "mcclust.ext", "fold")#,
             # .noexport = c("ldmvnorm_arma",
             #               "maketau",
             #               "makeSigma",
             #               "makeSigmacube",
             #               "makeSigmaDet",
             #               "mnorm_D_arma",
             #               "makeHellingerAvg")
             ) %dopar% {
                             clusterspirals(N=N)
                           }
sims <- matrix(unlist(x), ncol = 4, byrow=TRUE)
sims_means <- colMeans(sims)
round(sims_means,3)
sims_sd <- apply(sims, 2, sd)
round(sims_sd, 3)
