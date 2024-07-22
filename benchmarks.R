# libraries and preliminaries
require(fold)
require(mcclust)
require(mcclust.ext)
require(mclust)
require(dbscan)
require(dplyr)
require(readxl)
seed <- 1 # random number seed
### Example 1: Iris ###
## step 0: import data
c0 <- iris$Species
X <- cbind(iris$Sepal.Length, iris$Sepal.Width, iris$Petal.Length, iris$Petal.Width)
p <- ncol(X)
n <- nrow(X)

## step 1: compute FOLD
S <- 50000 # iterations
B <- 1000 # burnin
L <- 50
# fitting
set.seed(seed)
fit <- lmvnorm_gibbs(S = S,
                    y = scale(X),
                    L = L,
                    alpha = rep(1/2,L),
                    mu_0 = rep(0,p),
                    Sigma_0 = diag(1,p),
                    Sigma = 0.5*diag(1,p),
                    stops = 5000)
# burn-out
theta <- fit$theta[-(1:B)]
z <- fit$z[-(1:B),]
# thin
trip_ind <- seq(4,nrow(z),by=4)
theta <- theta[trip_ind]
z <- z[trip_ind,]
# computing Delta and candidates
Delta <- comp_delta(theta=theta, p=p, n = n)
cl <- hclust(as.dist(Delta), method = "average") # uses hierarchical clustering to get candidate set
# list of possible clusterings
max.k <- 10  # maximum number of clusters
candidates <- matrix(0, nrow = max.k, ncol = n)
for (h in 1:max.k){
  candidates[h,] = cutree(cl, k = h)
}
elbow(candidates,Delta)
c.fold <- candidates[2,]

## step 2: compute competitors
# W & G estimate -- VI loss
c.psm <- comp.psm(z)
c.mv <- minVI(psm = c.psm,
              cls.draw = z,
              max.k = max.k)
c.VI <- c.mv$cl
# W & G estimate -- Binder's Loss
c.b <- minbinder(psm = c.psm,
                 cls.draw = z,
                 max.k = max.k)
c.Binder <- c.b$cl
# Mclust
mcl <- Mclust(X, G = 1:max.k)
c.mcl <- mcl$classification
# K-means
km <- kmeans(X,centers=3)
c.km <- km$cluster
# HCA
hca <- hclust(dist(X,method="euclidean"), method="average")
c.hca <- cutree(hca, k=3)
# DBSCAN
kNNdistplot(X,k=ncol(X)) # suggests epsilon approx equal to 0.5
c.dbscan <- dbscan(X, minPts = 5, eps = 0.5)$cluster + 1 # for consistent labeling w/ other methods
# Spectral Clustering
D.iris = matrix(0, nrow=n, ncol = n)
bw <- 1 # bandwidth parameter
gauss.iris = exp(-(dist(X, method = "euclidean"))^2/bw^2)
D.iris[lower.tri(D.iris)] = gauss.iris
Dg = rowSums(D.iris + t(D.iris)) # Laplacian
Lg = diag(Dg) - (D.iris + t(D.iris))
evs = eigen(Lg) # eigenvalue decomposition
evals = rev(evs$values)
plot(evals[1:max.k]) # suggests 3 clusters
evecs = evs$vectors[,n:1]
k.spec <- 3
V <- evecs[,1:k.spec]
c.spec <- kmeans(V,centers=k.spec)$cluster

## step 3: save metrics
K_iris <- c(length(table(c.fold)),
             length(table(c.VI)),
             length(table(c.Binder)),
             length(table(c.mcl)),
             length(table(c.km)),
             length(table(c.hca)),
             length(table(c.dbscan)),
             length(table(c.spec)))
ARI_iris <- c(adjustedRandIndex(c0, c.fold),
               adjustedRandIndex(c0, c.VI),
               adjustedRandIndex(c0, c.Binder),
               adjustedRandIndex(c0, c.mcl),
               adjustedRandIndex(c0, c.km),
               adjustedRandIndex(c0, c.hca),
               adjustedRandIndex(c0, c.dbscan),
               adjustedRandIndex(c0, c.spec))


### Example 2: seeds ###

## step 0: import data
seeds <- read.table("data/seeds_dataset.txt")
c0 <- seeds$V8
X <- seeds %>% select(V3,V6,V7)
p <- ncol(X)
n <- nrow(X)

## step 1: compute FOLD
S <- 50000 # iterations
B <- 1000 # burnin
L <- 50
# fitting
set.seed(seed)
fit <- lmvnorm_gibbs(S = S,
                     y = scale(X),
                     L = L,
                     alpha = rep(1/2,L),
                     mu_0 = rep(0,p),
                     Sigma_0 = diag(1,p),
                     Sigma = 0.5*diag(1,p),
                     stops = 5000)
# burn-out
theta <- fit$theta[-(1:B)]
z <- fit$z[-(1:B),]
# thin
trip_ind <- seq(4,nrow(z),by=4)
theta <- theta[trip_ind]
z <- z[trip_ind,]
# computing Delta and candidates
Delta <- comp_delta(theta=theta, p=p, n = n)
cl <- hclust(as.dist(Delta), method = "average") # uses hierarchical clustering to get candidate set
# list of possible clusterings
max.k <- 10  # maximum number of clusters
candidates <- matrix(0, nrow = max.k, ncol = n)
for (h in 1:max.k){
  candidates[h,] = cutree(cl, k = h)
}
elbow(candidates,Delta)
c.fold <- candidates[3,]

## step 2: compute competitors
# W & G estimate -- VI loss
c.psm <- comp.psm(z)
c.mv <- minVI(psm = c.psm,
              cls.draw = z,
              max.k = max.k)
c.VI <- c.mv$cl
# W & G estimate -- Binder's Loss
c.b <- minbinder(psm = c.psm,
                 cls.draw = z,
                 max.k = max.k)
c.Binder <- c.b$cl
# Mclust
mcl <- Mclust(X, G = 1:max.k)
c.mcl <- mcl$classification
# K-means
km <- kmeans(X,centers=3)
c.km <- km$cluster
# HCA
hca <- hclust(dist(X,method="euclidean"), method="average")
c.hca <- cutree(hca, k=3)
# DBSCAN
kNNdistplot(X,k=ncol(X)) # suggests epsilon approx equal to 0.55
c.dbscan <- dbscan(X, minPts = 5, eps = 0.55)$cluster + 1 # for consistent labeling w/ other methods
# Spectral Clustering
D.seeds = matrix(0, nrow=n, ncol = n)
bw <- 1 # bandwidth parameter
gauss.seeds = exp(-(dist(X, method = "euclidean"))^2/bw^2)
D.seeds[lower.tri(D.seeds)] = gauss.seeds
Dg = rowSums(D.seeds + t(D.seeds)) # Laplacian
Lg = diag(Dg) - (D.seeds + t(D.seeds))
evs = eigen(Lg) # eigenvalue decomposition
evals = rev(evs$values)
plot(evals[1:max.k]) # suggests 2 clusters
evecs = evs$vectors[,n:1]
k.spec <- 2
V <- evecs[,1:k.spec]
c.spec <- kmeans(V,centers=k.spec)$cluster

## step 3: save metrics
K_seeds <- c(length(table(c.fold)),
             length(table(c.VI)),
             length(table(c.Binder)),
             length(table(c.mcl)),
             length(table(c.km)),
             length(table(c.hca)),
             length(table(c.dbscan)),
             length(table(c.spec)))
ARI_seeds <- c(adjustedRandIndex(c0, c.fold),
               adjustedRandIndex(c0, c.VI),
               adjustedRandIndex(c0, c.Binder),
               adjustedRandIndex(c0, c.mcl),
               adjustedRandIndex(c0, c.km),
               adjustedRandIndex(c0, c.hca),
               adjustedRandIndex(c0, c.dbscan),
               adjustedRandIndex(c0, c.spec))


### Example 3: Wine
## step 0: import data
wine <- read.csv("data/wine.data")
c0 <- wine$X1
Y <- wine[,-1]
X <- prcomp(Y, scale. = T)$x[,1:2]
p <- ncol(X)
n <- nrow(X)
## step 1: compute FOLD
S <- 50000 # iterations
B <- 1000 # burnin
L <- 50
# fitting
set.seed(seed)
fit <- lmvnorm_gibbs(S = S,
                     y = scale(X),
                     L = L,
                     alpha = rep(1/2,L),
                     mu_0 = rep(0,p),
                     Sigma_0 = diag(1,p),
                     Sigma = 0.5*diag(1,p),
                     stops = 5000)
# burn-out
theta <- fit$theta[-(1:B)]
z <- fit$z[-(1:B),]
# thin
trip_ind <- seq(4,nrow(z),by=4)
theta <- theta[trip_ind]
z <- z[trip_ind,]
# computing Delta and candidates
Delta <- comp_delta(theta=theta, p=p, n = n)
cl <- hclust(as.dist(Delta), method = "average") # uses hierarchical clustering to get candidate set
# list of possible clusterings
max.k <- 10  # maximum number of clusters
candidates <- matrix(0, nrow = max.k, ncol = n)
for (h in 1:max.k){
  candidates[h,] = cutree(cl, k = h)
}
elbow(candidates,Delta)
c.fold <- candidates[3,]

## step 2: compute competitors
# W & G estimate -- VI loss
c.psm <- comp.psm(z)
c.mv <- minVI(psm = c.psm,
              cls.draw = z,
              max.k = max.k)
c.VI <- c.mv$cl
# W & G estimate -- Binder's Loss
c.b <- minbinder(psm = c.psm,
                 cls.draw = z,
                 max.k = max.k)
c.Binder <- c.b$cl
# Mclust
mcl <- Mclust(X, G = 1:max.k)
c.mcl <- mcl$classification
# K-means
km <- kmeans(X,centers=3)
c.km <- km$cluster
# HCA
hca <- hclust(dist(X,method="euclidean"), method="average")
c.hca <- cutree(hca, k=3)
# DBSCAN
kNNdistplot(X,k=ncol(X)) # suggests epsilon approx equal to 0.7
c.dbscan <- dbscan(X, minPts = 5, eps = 0.7)$cluster + 1 # for consistent labeling w/ other methods
# Spectral Clustering
D.wine = matrix(0, nrow=n, ncol = n)
bw <- 1 # bandwidth parameter
gauss.wine = exp(-(dist(X, method = "euclidean"))^2/bw^2)
D.wine[lower.tri(D.wine)] = gauss.wine
Dg = rowSums(D.wine + t(D.wine)) # Laplacian
Lg = diag(Dg) - (D.wine + t(D.wine))
evs = eigen(Lg) # eigenvalue decomposition
evals = rev(evs$values)
plot(evals[1:max.k]) # suggests 4 clusters
evecs = evs$vectors[,n:1]
k.spec <- 4
V <- evecs[,1:k.spec]
c.spec <- kmeans(V,centers=k.spec)$cluster

## step 3: save metrics
K_wine <- c(length(table(c.fold)),
            length(table(c.VI)),
            length(table(c.Binder)),
            length(table(c.mcl)),
            length(table(c.km)),
            length(table(c.hca)),
            length(table(c.dbscan)),
            length(table(c.spec)))
ARI_wine <- c(adjustedRandIndex(c0, c.fold),
              adjustedRandIndex(c0, c.VI),
              adjustedRandIndex(c0, c.Binder),
              adjustedRandIndex(c0, c.mcl),
              adjustedRandIndex(c0, c.km),
              adjustedRandIndex(c0, c.hca),
              adjustedRandIndex(c0, c.dbscan),
              adjustedRandIndex(c0, c.spec))

### displaying metrics ###
K_metrics <- rbind(K_iris, K_seeds, K_wine)
ARI_metrics <- rbind(ARI_iris, ARI_seeds, ARI_wine)
colnames(K_metrics) <- colnames(ARI_metrics) <-
  c("FOLD", "VI", "Binder's", "EM", "K-Means", "HCA", "DBSCAN", "Spectral")
print(K_metrics)
print(round(ARI_metrics,3))

# ### Example 4: Breast Cancer Wisconsin
# ## step 0: import data
# wdbc <- read.csv("data/wdbc.data")
# c0 <- wdbc$M
# Y <- wdbc[,-(1:2)]
# X <- prcomp(Y, scale. = T)$x[,1:2]
# p <- ncol(X)
# n <- nrow(X)
# ## step 1: compute FOLD
# S <- 10000 # iterations
# B <- 500 # burnin
# L <- 50
# # fitting
# set.seed(seed)
# fit <- lmvnorm_gibbs(S = S,
#                      y = scale(X),
#                      L = L,
#                      alpha = rep(1/2,L),
#                      mu_0 = rep(0,p),
#                      Sigma_0 = diag(1,p),
#                      Sigma = 0.5*diag(1,p),
#                      stops = 5000)
# # burn-out
# theta <- fit$theta[-(1:B)]
# z <- fit$z[-(1:B),]
# # thin
# trip_ind <- seq(4,nrow(z),by=4)
# theta <- theta[trip_ind]
# z <- z[trip_ind,]
# # computing Delta and candidates
# Delta <- comp_delta(theta=theta, p=p, n = n)
# cl <- hclust(as.dist(Delta), method = "average") # uses hierarchical clustering to get candidate set
# # list of possible clusterings
# max.k <- 10  # maximum number of clusters
# candidates <- matrix(0, nrow = max.k, ncol = n)
# for (h in 1:max.k){
#   candidates[h,] = cutree(cl, k = h)
# }
# elbow(candidates,Delta)
# c.fold <- candidates[6,]
#
# ## step 2: compute competitors
# # W & G estimate -- VI loss
# c.psm <- comp.psm(z)
# c.mv <- minVI(psm = c.psm,
#               cls.draw = z,
#               max.k = max.k)
# c.VI <- c.mv$cl
# # W & G estimate -- Binder's Loss
# c.b <- minbinder(psm = c.psm,
#                  cls.draw = z,
#                  max.k = max.k)
# c.Binder <- c.b$cl
# # Mclust
# mcl <- Mclust(X, G = 1:max.k)
# c.mcl <- mcl$classification
# # K-means
# km <- kmeans(X,centers=3)
# c.km <- km$cluster
# # HCA
# hca <- hclust(dist(X,method="euclidean"), method="average")
# c.hca <- cutree(hca, k=3)
# # DBSCAN
# kNNdistplot(X,k=ncol(X)) # suggests epsilon approx equal to 0.4
# c.dbscan <- dbscan(X, minPts = 5, eps = 0.4)$cluster + 1 # for consistent labeling w/ other methods
# # Spectral Clustering
# D.wdbc = matrix(0, nrow=n, ncol = n)
# bw <- 1 # bandwidth parameter
# gauss.wdbc = exp(-(dist(X, method = "euclidean"))^2/bw^2)
# D.wdbc[lower.tri(D.wdbc)] = gauss.wdbc
# Dg = rowSums(D.wdbc + t(D.wdbc)) # Laplacian
# Lg = diag(Dg) - D.wdbc
# evs = eigen(Lg) # eigenvalue decomposition
# evals = rev(evs$values)
# plot(evals[1:max.k]) # suggests 6 clusters
# evecs = evs$vectors[,n:1]
# k.spec <- 6
# V <- evecs[,1:k.spec]
# c.spec <- kmeans(V,centers=k.spec)$cluster
#
# ## step 3: save metrics
# K_wdbc <- c(length(table(c.fold)),
#             length(table(c.VI)),
#             length(table(c.Binder)),
#             length(table(c.mcl)),
#             length(table(c.km)),
#             length(table(c.hca)),
#             length(table(c.dbscan)),
#             length(table(c.spec)))
# ARI_wdbc <- c(adjustedRandIndex(c0, c.fold),
#               adjustedRandIndex(c0, c.VI),
#               adjustedRandIndex(c0, c.Binder),
#               adjustedRandIndex(c0, c.mcl),
#               adjustedRandIndex(c0, c.km),
#               adjustedRandIndex(c0, c.hca),
#               adjustedRandIndex(c0, c.dbscan),
#               adjustedRandIndex(c0, c.spec))
#
# ### displaying metrics ###
# K_metrics <- rbind(K_iris, K_seeds, K_wine, K_wdbc)
# ARI_metrics <- rbind(ARI_iris, ARI_seeds, ARI_wine, K_wdbc)
#
