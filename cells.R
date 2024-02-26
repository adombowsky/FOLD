### file to carry out analysis on GSE81861 Cell Line Dataset ###
# libraries
require(sn)
require(ggplot2)
require(mcclust)
require(mcclust.ext)
require(mclust)
require(grid)
require(gridExtra)
require(Rcpp)
require(RcppArmadillo)
require(plyr)
require(dplyr)
require(reshape2)
require(mvtnorm)
require(scran)
require(scuttle)
require(M3Drop)
require(latex2exp)
require(cowplot)
require(uwot)
source("rfuncts/simple_SALSO.R")
source("rfuncts/salso_parallel.R")
sourceCpp("rcppfuncts/salso.cpp")
set.seed(2001)
# loading data and preprocessing
cells_csv <- read.csv("data/cells.csv")
cells <- as.matrix(cells_csv[,-1])
## remove dimnames
cells <- SingleCellExperiment(list(counts=cells))
## remove low count cells
qcstats <- perCellQCMetrics(cells)
qcfilter <- quickPerCellQC(qcstats) # come back to
cells <- cells[,!qcfilter$discard]
## perform SCRAN normalization
clusters <- quickCluster(cells)
cells <- computeSumFactors(cells, clusters=clusters)
cells <- logNormCounts(cells, log = T) 
data <- cells@assays@data@listData$logcounts
norm_data <- M3DropConvertData(data,
                               is.log = T)
norm_data <- as.data.frame(norm_data)
cellnames <- colnames(norm_data)
norm_data <- unname(norm_data)
drop_data <- M3DropFeatureSelection(norm_data, mt_method="fdr", mt_threshold=0.01)
genes <- drop_data$Gene
data <- norm_data[genes,]
n <- ncol(data)
## configuring names
dictionary <- c("_A549_","_GM12878_", "_H1437_", "_HCT116_", "_IMR90_", "_H1_", "_K562_")
K <- length(dictionary)
types <- matrix(0, nrow = n, ncol = K)
for (k in 1:K) {
  types[,k] <- k * grepl(dictionary[k], cellnames)
}
types <- rowSums(types)
## taking PCA
d <- 5
dpc <- prcomp(t(data), rank. = d, scale. = TRUE)
X <- scale(dpc$x)

# running MCMC
source("rfuncts/mvnorm_gibbs.R")
source("rfuncts/mnorm_D_apply.R")
source("rfuncts/foldball.R")
sourceCpp("rcppfuncts/mnorm_D_arma.cpp")
S <- 25000 # iterations
B <- 1000 # burnin
L <- 50
# fitting
fit <- mvnorm_gibbs(S = S,
                    y = X,
                    L = L,
                    alpha = rep(1/2,L),
                    w = rep(0,d),
                    kappa = 1,
                    r = d+2,
                    C = diag(1,d),
                    stops = 1000)
# saveRDS(fit, file = "output/cells/rna_fit.Rdata")
#fit <- readRDS("output/cells/rna_fit.Rdata")
# burn-out
theta <- fit$theta[-(1:B)]
M <- length(theta)
n <- nrow(X)
theta <- array(unlist(theta), dim = c(n, 2*d + choose(d,2), M))
z <- fit$z[-(1:B),]
# thin
trip_ind <- seq(4,M,by=4)
theta <- theta[,,trip_ind]
z <- z[trip_ind,]

# apply hellinger distance function
Delta_UT <- makeHellingerAvg(theta = theta, d = d, n = n)
#write.csv(Delta_UT, "output/delta_ut_cells.csv", row.names = F)
#Delta_UT <- read.csv("output/delta_ut_cells.csv")$V1
Delta <- matrix(0, nrow = n, ncol = n)
Delta[upper.tri(Delta)] <- Delta_UT
Delta <- Delta + t(Delta)
Delta_H <- Delta

# choose candidate set
source("rfuncts/mergeloss.R")
# single-linkage
Delta <- as.dist(Delta)
cl <- hclust(Delta, method = "average") # uses hierarchical clustering to get candidate set
# list of possible clusterings
max.k <- 24  # maximum number of clusters
labs <- matrix(0, nrow = max.k, ncol = n)
for (h in 1:max.k){
  labs[h,] = cutree(cl, k = h)
}
w_avg <- sum(Delta_UT)/(sum(1 - Delta_UT))
gamma_med <- median(Delta_UT)
w_med = gamma_med/(1-gamma_med)
losses <- apply(labs, 1, mergeloss, Delta_lt = Delta_H[lower.tri(Delta_H)], w=w_avg)
# minimize
min_loss <- which.min(losses)
c.fold = labs[min_loss,]
# W & G estimate -- VI loss
c.psm <- comp.psm(z)
c.mv <- minVI(psm = c.psm,
              cls.draw = z,
              max.k = max.k)
c.VI <- c.mv$cl
#write.csv(c.VI, "output/cVI_cells.csv", row.names = F)
# W & G estimate -- Binder's Loss
c.b <- minbinder(psm = c.psm,
                 cls.draw = z,
                 max.k = max.k)
c.Binder <- c.b$cl
#write.csv(c.Binder, "output/cBinder_cells.csv", row.names = F)
# Elbow Plot
W <- seq(0,250,by=2)
tv <- c()
tv <- c()
h_bar <- sum(Delta_H[lower.tri(Delta_H)])
for (k in 1:max.k){
  wss <- c()
  for (h in 1:k) {
    Delta_h <- Delta_H[labs[k,]==h, labs[k,]==h]
    wss[h] <- sum(Delta_h[lower.tri(Delta_h)])
  }
  tv[k] <- sum(wss)/h_bar
}
elb_df <- data.frame(Size = 1:max.k,
                     TV = tv)
ggplot(elb_df, aes(x = Size, y  = TV)) + geom_line() + geom_point() + 
  labs(title = "Elbow Plot") +
  xlab("No. of Clusters") +
  ylab(latex2exp::TeX("$r_{\\omega}$")) +
  scale_x_discrete(limits = as.character(1:max.k)) +
  theme_bw() +
  theme(text = element_text(size=18))
c.fold = labs[6,]
# alternatively, using SALSO
c.salso = salso_parallel(M=150,
                         Delta=Delta_H,
                         omega = 25)
## compare risk value
mergeloss(c.fold, Delta_lt = Delta_H[lower.tri(Delta_H)], w=25)/choose(n,2)
mergeloss(c.salso, Delta_lt = Delta_H[lower.tri(Delta_H)], w=25)/choose(n,2) # salso generally gets higher risk value
# Mclust
mcl <- Mclust(X, G = 1:24)
# plot with UMAP
set.seed(2001)
umapped_cells <- uwot::umap(t(data), 
                            min_dist = 0.8)
types <- as.factor(types)
levels(types) <- c("A5A9", "GM12878", "H1437", "HCT116",
                   "IMR90", "H1", "K562")
umapped_df <- data.frame(umapped_cells,
                         c.VI = factor(c.VI),
                         c.Binder <- factor(c.Binder),
                         c.Mclust= factor(mcl$classification),
                         c.fold = factor(c.fold),
                         types = factor(types))
colnames(umapped_df)[1:2] <- c("UMAP1", "UMAP2")
gg_truth <- umapped_df %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = types)) + geom_point() + 
  labs(caption = "The true cell types.", color = "Types", title = " ") + theme_bw() +
  theme(
    #axis.line = element_blank(),
    #    axis.ticks = element_blank(),
    #    axis.text = element_blank(),
    text = element_text(size=14),
    plot.background = element_blank(),
    plot.caption = element_text(hjust=0.5)
    #    panel.grid.major = element_blank(),
    #    panel.grid.minor = element_blank(),
    #    panel.border = element_blank()
  )
gg_VI <- umapped_df %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = c.VI)) + geom_point() + 
  labs(caption = "Clusters given by minimizing the VI and Binder's loss.", color = "Cluster", title = " ") + theme_bw() +
  theme(
    #axis.line = element_blank(),
    #    axis.ticks = element_blank(),
    #    axis.text = element_blank(),
    text = element_text(size=14),
    plot.background = element_blank(),
    plot.caption = element_text(hjust=0.5)
    #    panel.grid.major = element_blank(),
    #    panel.grid.minor = element_blank(),
    #    panel.border = element_blank()
  )
gg_mclust <- umapped_df %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = c.Mclust)) + geom_point() + 
  labs(caption = "Clusters given by Mclust.", color = "Cluster", title = " ") + theme_bw() +
  theme(
    #axis.line = element_blank(),
    #    axis.ticks = element_blank(),
    #    axis.text = element_blank(),
    text = element_text(size=14),
    plot.background = element_blank(),
    plot.caption = element_text(hjust=0.5)
    #    panel.grid.major = element_blank(),
    #    panel.grid.minor = element_blank(),
    #    panel.border = element_blank()
  )
gg_fold <- umapped_df %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = c.fold)) + geom_point() + 
  labs(caption = latex2exp::TeX("Clusters given by FOLD."), color = "Cluster", title = " ") + theme_bw() +
  theme(
    #axis.line = element_blank(),
    #    axis.ticks = element_blank(),
    #    axis.text = element_blank(),
    text = element_text(size=14),
    plot.background = element_blank(),
    plot.caption = element_text(hjust=0.5)
    #    panel.grid.major = element_blank(),
    #    panel.grid.minor = element_blank(),
    #    panel.border = element_blank()
  )
plot_grid(gg_truth, gg_fold, gg_VI, gg_mclust, labels = c('(a)', '(b)', '(c)', '(d)'))
# hellinger distance matrix
## re-label c.fold by size
library(reshape2)
library(latex2exp)
c.fold.ordered <- order(c.fold)
Delta.ordered <- 1-Delta_H[c.fold.ordered, c.fold.ordered]
melt.delta <- melt(Delta.ordered)
hdmat = ggplot(data = melt.delta, aes(x = Var1, y = rev(Var2), fill = value)) +
  geom_tile() + scale_fill_distiller(palette = "YlOrBr")  +
  xlab("") + ylab("") + labs(title = "Hellinger Similarity Matrix",
                             fill = TeX("$1-\\Delta_{ij}$")) +
  theme_bw() + 
  theme(axis.line = element_blank(),
        legend.position = "left",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
# posterior similarity matrix
c.VI.ordered <- order(c.VI)
psm.ordered <- c.psm[c.VI.ordered, c.VI.ordered]
melt.psm <- melt(psm.ordered)
psm = ggplot(data = melt.psm, aes(x = Var1, y = rev(Var2), fill = value)) +
  geom_tile() + scale_fill_distiller(palette = "YlOrBr") +
  xlab("") + ylab("") + labs(title = "Posterior Similiarity Matrix",
                             fill = TeX("$\\Pi(s_i = s_j | X)$")) +
  theme_bw() + 
  theme(axis.line = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
grid.arrange(hdmat, psm, nrow = 1)
# credible ball and cGPSM
#cgsamps <- getcGsamps(theta=theta, w=25, d=5, max.k=10)
cgsamps <- read.csv("output/cells/cgsamps.csv")
cgsamps <- as.matrix(cgsamps)
cb <- cGball(c.fold,cgsamps)
foldpsm <- cGpsm(cgsamps)
# plot PSM
c.fold.ordered <- order(c.fold)
fpsm.ordered <- foldpsm[c.fold.ordered, c.fold.ordered]
melt.psm <- melt(fpsm.ordered)
fpsm = ggplot(data = melt.psm, aes(x = Var1, y = rev(Var2), fill = value)) +
  geom_tile() + scale_fill_distiller(palette = "YlOrBr") +
  xlab("") + ylab("") + labs(title = "Posterior Similiarity Matrix",
                             fill = TeX("$\\Pi(s_i = s_j | X)$")) +
  theme_bw() + 
  theme(axis.line = element_blank(),
        legend.position = "right",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
# credible ball partitions
ball.df <- data.frame(
  umapped_cells,
  cfold = factor(c.fold),
  choriz = factor(cb$c.horiz[1,]),
  c.upperv = factor(cb$c.uppervert[1,]),
  c.lowerv = factor(cb$c.lowervert[1,]))
foldumap <- ggplot(ball.df, aes(x = X1, y = X2, color = cfold)) + geom_point() + labs(color = "Cluster", caption = latex2exp::TeX("The point estimator, $\\textbf{c}_{\\textrm{FOLD}}.$"), 
                                                                                      title = " ")  +
  xlab("UMAP1") + ylab("UMAP2") +
  theme_bw() + 
  #scale_color_brewer(palette="Dark2")+
  theme(
    #axis.line = element_blank(),
    #    axis.ticks = element_blank(),
    #    axis.text = element_blank(),
    text = element_text(size=14),
    plot.background = element_blank(),
    plot.caption = element_text(hjust = 0.5)
    #    panel.grid.major = element_blank(),
    #    panel.grid.minor = element_blank(),
    #    panel.border = element_blank()
  ) 
horizumap <- ggplot(ball.df, aes(x = X1, y = X2, color = choriz)) + geom_point() + labs(caption="95% credible ball horizontal bound.", color = "Cluster", 
                                                                                        title = " ") +
  xlab("UMAP1") + ylab("UMAP2") +
  theme_bw() + 
  #scale_color_brewer(palette="Dark2")+
  theme(
    #axis.line = element_blank(),
    #    axis.ticks = element_blank(),
    #    axis.text = element_blank(),
    text = element_text(size=14),
    plot.background = element_blank(),
    plot.caption = element_text(hjust = 0.5)
    #    panel.grid.major = element_blank(),
    #    panel.grid.minor = element_blank(),
    #    panel.border = element_blank()
  )
uppervumap <- ggplot(ball.df, aes(x = X1, y = X2, color = c.upperv)) + geom_point() + labs(caption="95% credible ball vertical upper bound.", color = "Cluster",
                                                                                           title = " ")  +
  xlab("UMAP1") + ylab("UMAP2") +
  theme_bw() + 
  #scale_color_brewer(palette="Dark2")+
  theme(
    #axis.line = element_blank(),
    #    axis.ticks = element_blank(),
    #    axis.text = element_blank(),
    text = element_text(size=14),
    plot.background = element_blank(),
    plot.caption = element_text(hjust = 0.5)
    #    panel.grid.major = element_blank(),
    #    panel.grid.minor = element_blank(),
    #    panel.border = element_blank()
  )
lowervumap <- ggplot(ball.df, aes(x = X1, y = X2, color = c.lowerv)) + geom_point() + labs(caption="95% credible ball vertical lower bound.", color = "Cluster", title = " ")  +
  xlab("UMAP1") + ylab("UMAP2") +
  theme_bw() + 
  #scale_color_brewer(palette="Dark2")+
  theme(
    #axis.line = element_blank(),
    #    axis.ticks = element_blank(),
    #    axis.text = element_blank(),
    text = element_text(size=14),
    plot.background = element_blank(),
    plot.caption = element_text(hjust = 0.5)
    #    panel.grid.major = element_blank(),
    #    panel.grid.minor = element_blank(),
    #    panel.border = element_blank()
  )
plot_grid(foldumap, horizumap, lowervumap, uppervumap, labels = c('(a)', '(b)', '(c)', '(d)'))

# algorithmic methods
set.seed(2001)
## competitor: hierarchical clustering
dist.cells = dist(X, method = "euclidean")
h.cells = hclust(dist.cells, method = "average")
c.hc = cutree(h.cells,k=7) # with the correctly specified number of clusters
umapped_df$c.hc <- as.factor(c.hc)
gg_hc <- umapped_df %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = c.hc)) + geom_point() + 
  labs(caption = "Avg. Linkage Hierarchical Clustering.", color = "Cluster", title = " ") + theme_bw() +
  theme(
    #axis.line = element_blank(),
    #    axis.ticks = element_blank(),
    #    axis.text = element_blank(),
    text = element_text(size=14),
    plot.background = element_blank(),
    plot.caption = element_text(hjust=0.5)
    #    panel.grid.major = element_blank(),
    #    panel.grid.minor = element_blank(),
    #    panel.border = element_blank()
  )

## competitor: k-means
### first, make elbow plot
k.max <- 10
tws <- c()
for (h in 1:k.max){
  km <- kmeans(x=X, centers = h)
  tws[h] <- km$tot.withinss
}
# Generate a data frame containing both k and tot_withinss
elbow_df <- data.frame(
  k = 1:k.max,
  tot_withinss = tws
)
# Plot the elbow plot
ggplot(elbow_df, aes(x = k, y = tot_withinss)) +
  geom_line() + geom_point()+
  scale_x_continuous(breaks = 1:10)
km.cells = kmeans(X,centers=7)
c.km = km.cells$cluster
umapped_df$c.km <- as.factor(c.km)
gg_km <- umapped_df %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = c.km)) + geom_point() + 
  labs(caption = "K-means.", color = "Cluster", title = " ") + theme_bw() +
  theme(
    #axis.line = element_blank(),
    #    axis.ticks = element_blank(),
    #    axis.text = element_blank(),
    text = element_text(size=14),
    plot.background = element_blank(),
    plot.caption = element_text(hjust=0.5)
    #    panel.grid.major = element_blank(),
    #    panel.grid.minor = element_blank(),
    #    panel.border = element_blank()
  )
# competitor: dbscan
library(dbscan)
# selecting epsilon
kNNdistplot(X,k=ncol(X)) # suggests epsilon approx equal to 1
cells.dbs <- dbscan(X, minPts = ncol(X)+1,eps=1)
c.dbs = cells.dbs$cluster + 1 # for consistent labeling w/ other methods
umapped_df$c.dbs = factor(c.dbs)
gg_dbs <- umapped_df %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = c.dbs)) + geom_point() + 
  labs(caption = "DBSCAN.", color = "Cluster", title = " ") + theme_bw() +
  theme(
    #axis.line = element_blank(),
    #    axis.ticks = element_blank(),
    #    axis.text = element_blank(),
    text = element_text(size=14),
    plot.background = element_blank(),
    plot.caption = element_text(hjust=0.5)
    #    panel.grid.major = element_blank(),
    #    panel.grid.minor = element_blank(),
    #    panel.border = element_blank()
  )

## competitor: spectral clustering
k.spec.max = 10
D.cells = matrix(0, nrow=n, ncol = n)
dist.cells = dist(X, method = "euclidean")
bw <- 1 # bandwidth parameter
gauss.cells = exp(-(dist.cells)^2/bw^2)
D.cells[lower.tri(D.cells)] = gauss.cells
D.cells = D.cells + t(D.cells)
# Laplacian 
Dg = rowSums(D.cells)
Lg = diag(Dg) - D.cells
# eigenvalue decomposition
evs = eigen(Lg)
evals = rev(evs$values)
plot(evals[1:k.spec.max])
evecs = evs$vectors[,n:1]
k.spec <- 5
V <- evecs[,1:k.spec]

# now do k-means
km.sc <- kmeans(V,centers=k.spec)
c.sc <- km.sc$cluster
umapped_df$c.sc = factor(c.sc)
gg_sc <- umapped_df %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = c.sc)) + geom_point() + 
  labs(caption = "Spectral Clustering.", color = "Cluster", title = " ") + theme_bw() +
  theme(
    #axis.line = element_blank(),
    #    axis.ticks = element_blank(),
    #    axis.text = element_blank(),
    text = element_text(size=14),
    plot.background = element_blank(),
    plot.caption = element_text(hjust=0.5)
    #    panel.grid.major = element_blank(),
    #    panel.grid.minor = element_blank(),
    #    panel.border = element_blank()
  )

# comparing all methods + plots for algorithmic methods
# gives number of clusters and adjusted Rand w/ truth
summary.df <- data.frame(FOLD=c(adjustedRandIndex(c.fold,types), length(table(c.fold))),
                         VI = c(adjustedRandIndex(c.VI,types), length(table(c.VI))),
                         Mclust = c(adjustedRandIndex(factor(mcl$classification),types), length(table(factor(mcl$classification)))),
                         HC = c(adjustedRandIndex(c.hc,types), length(table(c.hc))),
                         kmeans = c(adjustedRandIndex(c.km,types), length(table(c.km))),
                         dbscan = c(adjustedRandIndex(c.dbs,types), length(table(c.dbs))),
                         spectral = c(adjustedRandIndex(c.sc,types), length(table(c.sc)))
)
xtable::xtable(summary.df,digits = 3)
# plotting all algorithmic methods together (for the supplement)
plot_grid(gg_hc, gg_km, gg_dbs, gg_sc, labels =c("(a)", "(b)", "(c)", "(d)"))
