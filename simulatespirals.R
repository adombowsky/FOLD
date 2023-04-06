### simulation study given in the Supplementary Material ###
# see the file "clusterspirals.R" in the "rfuncts" folder for one replication 
# packages
library(grid)
library(gridExtra)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doParallel)
library(plyr)
library(dplyr)
library(reshape2)
# r functions
source("rfuncts/mvnorm_gibbs.R")
source("rfuncts/mnorm_D_apply.R")
source("rfuncts/mergeloss.R")
source("rfuncts/clusterspirals.R")
source("rfuncts/foldball.R")
sourceCpp("rcppfuncts/mnorm_D_arma.cpp")
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
             .packages = c("mcclust", "Rcpp", "cluster", "KODAMA", "mcclust.ext"),
             .noexport = c("ldmvnorm_arma",
                           "maketau",
                           "makeSigma",
                           "makeSigmacube",
                           "makeSigmaDet",
                           "mnorm_D_arma",
                           "makeHellingerAvg")) %dopar% {
                             clusterspirals(N=N)
                           }
sims <- matrix(unlist(x), ncol = 3, byrow=TRUE)
sims_means <- colMeans(sims)
round(sims_means,3)
sims_sd <- apply(sims, 2, sd)
round(sims_sd, 3)
