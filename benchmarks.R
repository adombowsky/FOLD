# libraries and preliminaries
require(foldcluster)
require(mcclust)
require(mcclust.ext)
require(mclust)
require(dbscan)
require(dplyr)
require(readxl)
require(ggplot2)
require(cowplot)
seed <- 1 # random number seed
### Example 1: Iris ###
## step 0: import data
c0 <- iris$Species
c0 <- 1*(c0=="setosa") + 2*(c0!="setosa") # merging versicolor and setosa
X <- cbind(iris$Sepal.Length, iris$Sepal.Width, iris$Petal.Length, iris$Petal.Width)
p <- ncol(X)
n <- nrow(X)

## step 1: compute FOLD
S <- 50000 # iterations
B <- 1000 # burnin
L <- 50
# fitting
set.seed(seed)
fit <- mvnorm_gibbs(S = S,
                    y = scale(X),
                    L = L,
                    alpha = rep(1/2,L),
                    w=rep(0,p),
                    kappa=1,
                    r=p+2,
                    C=diag(1,p),
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
km <- kmeans(X,centers=2)
c.km <- km$cluster
# HCA
hca <- hclust(dist(X,method="euclidean"), method="average")
c.hca <- cutree(hca, k=2)
# DBSCAN
kNNdistplot(X,k=ncol(X)) # suggests epsilon approx equal to 0.6
c.dbscan <- dbscan(X, minPts = 5, eps = 0.6)$cluster + 1 # for consistent labeling w/ other methods
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
k.spec <- 2
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

### Example 2: Flea-Beetles

## step 0: import data
# Dataset link: http://www.dm.unibo.it/~simoncin/FleaBeetles.html
# move into "data" folder
fleabeetles <- read.delim("data/flea_beetles")
c0 <- as.factor(fleabeetles$Species)
X <- fleabeetles[,-3]
p <- ncol(X)
n <- nrow(X)
## step 1: compute FOLD
S <- 50000 # iterations
B <- 1000 # burnin
L <- 50
# fitting
set.seed(seed)
fit <- mvnorm_gibbs(S = S,
                     y = scale(X),
                     L = L,
                     alpha = rep(1/2,L),
                     w=rep(0,p),
                     kappa=1,
                     r=p+2,
                     C=diag(1,p),
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
c.dbscan <- dbscan(X, minPts = 5, eps = 1)$cluster + 1 # for consistent labeling w/ other methods
# Spectral Clustering
D.fb = matrix(0, nrow=n, ncol = n)
bw <- 1 # bandwidth parameter
gauss.fb = exp(-(dist(X, method = "euclidean"))^2/bw^2)
D.fb[lower.tri(D.fb)] = gauss.fb
Dg = rowSums(D.fb + t(D.fb)) # Laplacian
Lg = diag(Dg) - (D.fb + t(D.fb))
evs = eigen(Lg) # eigenvalue decomposition
evals = rev(evs$values)
plot(evals[1:max.k]) # suggests 4 clusters
evecs = evs$vectors[,n:1]
k.spec <- 5
V <- evecs[,1:k.spec]
c.spec <- kmeans(V,centers=k.spec)$cluster

## step 3: save metrics
K_fb <- c(length(table(c.fold)),
          length(table(c.VI)),
          length(table(c.Binder)),
          length(table(c.mcl)),
          length(table(c.km)),
          length(table(c.hca)),
          length(table(c.dbscan)),
          length(table(c.spec)))
ARI_fb <- c(adjustedRandIndex(c0, c.fold),
            adjustedRandIndex(c0, c.VI),
            adjustedRandIndex(c0, c.Binder),
            adjustedRandIndex(c0, c.mcl),
            adjustedRandIndex(c0, c.km),
            adjustedRandIndex(c0, c.hca),
            adjustedRandIndex(c0, c.dbscan),
            adjustedRandIndex(c0, c.spec))


### Example 3: Wine

## step 0: import data
# Dataset link: https://archive.ics.uci.edu/dataset/109/wine
# move into "data" folder
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
fit <- mvnorm_gibbs(S = S,
                    y = scale(X),
                    L = L,
                    alpha = rep(1/2,L),
                    w=rep(0,p),
                    kappa=1,
                    r=p+2,
                    C=diag(1,p),
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
K_metrics <- rbind(K_iris, K_fb, K_wine)
ARI_metrics <- round(rbind(ARI_iris, ARI_fb, ARI_wine),3)
colnames(K_metrics) <- colnames(ARI_metrics) <-
  c("FOLD", "VI", "Binder's", "EM", "K-Means", "HCA", "DBSCAN", "Spectral")
print(K_metrics)
print(ARI_metrics)


###
### Univariate Examples ###
###

### Example 4: Galaxy

## step 0: import data
# Dataset link: https://people.maths.bris.ac.uk/~mapjg/mixdata
# move into "data" folder
X <- read.delim("data/galaxy")$X82
p <- 1
n <- length(X)
## step 1: compute FOLD
S <- 50000 # iterations
B <- 1000 # burnin
L <- 10
# fitting
set.seed(seed)
fit <- unorm_gibbs(S=S,
                   y=as.vector(scale(X)),
                   L=L,
                   alpha = rep(1/2,L),
                   w = 0,
                   kappa = 1,
                   r = 2,
                   C = 2,
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

# plotting
galaxy_df <- data.frame(data = X,
                        FOLD = as.factor(c.fold),
                        VI = as.factor(c.VI),
                        Binder = as.factor(c.Binder),
                        EM = as.factor(c.mcl))
# ggplot(galaxy_df, aes(x=X)) +
#   geom_histogram() + theme_bw() # histogram of whole dataset
galaxy.fold <- ggplot(galaxy_df, aes(x=X, fill = FOLD)) +
  geom_histogram() + theme_bw() +
  ylab("Count") + xlab("Data") +
  theme(text = element_text(size=15))
galaxy.vi <- ggplot(galaxy_df, aes(x=X, fill = VI)) +
  geom_histogram() + theme_bw() +
  ylab("Count") + xlab("Data") +
  theme(text = element_text(size=15))
galaxy.binder <- ggplot(galaxy_df, aes(x=X, fill = Binder)) +
  geom_histogram() + theme_bw() +
  ylab("Count") + xlab("Data") +
  theme(text = element_text(size=15))
galaxy.em <- ggplot(galaxy_df, aes(x=X, fill = EM)) +
  geom_histogram() + theme_bw() +
  ylab("Count") + xlab("Data") +
  theme(text = element_text(size=15))
galaxy_grid <- plot_grid(galaxy.fold, galaxy.vi, galaxy.binder, galaxy.em, nrow = 2,
          labels = c("(a)", "(b)", "(c)", "(d)"))
galaxy_title <- ggdraw() +
  draw_label(
    "Galaxy Data",
    fontface = 'bold',
    x = 0,
    hjust = 0,
    size = 20
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 208)
  )
plot_grid(
  galaxy_title, galaxy_grid,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.07, 1)
)

### Example 5: Acidity

## step 0: import data
# Dataset link: https://people.maths.bris.ac.uk/~mapjg/mixdata
# move into "data" folder
X <- read.delim("data/acidity")$X155
p <- 1
n <- length(X)
## step 1: compute FOLD
S <- 50000 # iterations
B <- 1000 # burnin
L <- 10
# fitting
set.seed(seed)
fit <- unorm_gibbs(S=S,
                   y=as.vector(scale(X)),
                   L=L,
                   alpha = rep(1/2,L),
                   w = 0,
                   kappa = 1,
                   r = 2,
                   C = 2,
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

# plotting
acidity_df <- data.frame(data = X,
                        FOLD = as.factor(c.fold),
                        VI = as.factor(c.VI),
                        Binder = as.factor(c.Binder),
                        EM = as.factor(c.mcl))
acidity.fold <- ggplot(acidity_df, aes(x=X, fill = FOLD)) +
  geom_histogram() + theme_bw() +
  ylab("Count") + xlab("Data") +
  theme(text = element_text(size=15))
acidity.vi <- ggplot(acidity_df, aes(x=X, fill = VI)) +
  geom_histogram() + theme_bw() +
  ylab("Count") + xlab("Data") +
  theme(text = element_text(size=15))
acidity.binder <- ggplot(acidity_df, aes(x=X, fill = Binder)) +
  geom_histogram() + theme_bw() +
  ylab("Count") + xlab("Data") +
  theme(text = element_text(size=15))
acidity.em <- ggplot(acidity_df, aes(x=X, fill = EM)) +
  geom_histogram() + theme_bw() +
  ylab("Count") + xlab("Data") +
  theme(text = element_text(size=15))
acidity_grid <- plot_grid(acidity.fold, acidity.vi, acidity.binder, acidity.em, nrow = 2,
          labels = c("(a)", "(b)", "(c)", "(d)"))
acidity_title <- ggdraw() +
  draw_label(
    "Acidity Data",
    fontface = 'bold',
    x = 0,
    hjust = 0,
    size = 20
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 208)
  )
plot_grid(
  acidity_title, acidity_grid,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.07, 1)
)

### Example 6: Enzyme

## step 0: import data
# Dataset link: https://people.maths.bris.ac.uk/~mapjg/mixdata
# move into "data" folder
X <- read.delim("data/enzyme")$X245
p <- 1
n <- length(X)
## step 1: compute FOLD
S <- 50000 # iterations
B <- 1000 # burnin
L <- 10
# fitting
set.seed(seed)
fit <- unorm_gibbs(S=S,
                   y=as.vector(scale(X)),
                   L=L,
                   alpha = rep(1/2,L),
                   w = 0,
                   kappa = 1,
                   r = 2,
                   C = 2,
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

# plotting
enzyme_df <- data.frame(data = X,
                         FOLD = as.factor(c.fold),
                         VI = as.factor(c.VI),
                         Binder = as.factor(c.Binder),
                         EM = as.factor(c.mcl))
enzyme.fold <- ggplot(enzyme_df, aes(x=X, fill = FOLD)) +
  geom_histogram() + theme_bw() +
  ylab("Count") + xlab("Data") +
  theme(text = element_text(size=15))
enzyme.vi <- ggplot(enzyme_df, aes(x=X, fill = VI)) +
  geom_histogram() + theme_bw() +
  ylab("Count") + xlab("Data") +
  theme(text = element_text(size=15))
enzyme.binder <- ggplot(enzyme_df, aes(x=X, fill = Binder)) +
  geom_histogram() + theme_bw() +
  ylab("Count") + xlab("Data") +
  theme(text = element_text(size=15))
enzyme.em <- ggplot(enzyme_df, aes(x=X, fill = EM)) +
  geom_histogram() + theme_bw() +
  ylab("Count") + xlab("Data") +
  theme(text = element_text(size=15))
enzyme_grid <- plot_grid(enzyme.fold, enzyme.vi, enzyme.binder, enzyme.em, nrow = 2,
          labels = c("(a)", "(b)", "(c)", "(d)"))
enzyme_title <- ggdraw() +
  draw_label(
    "Enzyme Data",
    fontface = 'bold',
    x = 0,
    hjust = 0,
    size = 20
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 208)
  )
plot_grid(
  enzyme_title, enzyme_grid,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(0.07, 1)
)

