### file to make all figures in main paper related to moons data ###
# packages
library(sn)
library(ggplot2)
library(dplyr)
library(mcclust.ext)
library(RColorBrewer)
library(ggsci)
library(RSSL)
library(reshape2)
library(ggforce)
library(Rcpp)
library(RcppArmadillo)
library(cowplot)
# moons data
set.seed(1999)
d <- 2
data <- generateCrescentMoon(n=500, d=d, sigma=0.5)
y <- scale(as.matrix(data[,-1]))
# fitting a GMM
source("rfuncts/locationnorm_gibbs.R")
source("rfuncts/mnorm_D_apply.R")
source("rfuncts/foldball.R")
source("rfuncts/loc_foldball.R")
sourceCpp("rcppfuncts/locnorm_D_arma.cpp")
S <- 35000 # iterations
B <- 1000 # burnin
L <- 30 # components
sig_sq <- 0.02 # diagonal covariance term
Sigma <- sig_sq*diag(1,d)
# fitting
fit <- loc_mvnorm_gibbs(S = S,
                    y = y,
                    L = L,
                    alpha = rep(1/L,L),
                    mu_0 = rep(0,d),
                    Sigma_0 = diag(2,d),
                    Sigma = Sigma,
                    stops = 50)
# burn-out
theta <- fit$theta[-(1:B)]
M <- length(theta)
n <- nrow(y)
theta <- array(unlist(theta), dim = c(n, d, M))
z <- fit$z[-(1:B),]
# thin
trip_ind <- seq(4,M,by=4)
theta <- theta[,,trip_ind]
z <- z[trip_ind,]
# cluster
# VI
c.psm <- mcclust::comp.psm(z)
c.mv <- mcclust.ext::minVI(psm = c.psm,
                           cls.draw = z,
                           max.k = 10)
c.VI <- c.mv$cl
# Binder's
c.mv <- mcclust::minbinder(psm = c.psm,
                           cls.draw = z,
                           max.k = 10)
c.B <- c.mv$cl
# FOLD
Delta_UT <- locHellingerAvg(theta = theta, Sigma = Sigma, d = d, n = n)
Delta <- matrix(0, nrow = n, ncol = n)
Delta[upper.tri(Delta)] <- Delta_UT
Delta <- Delta + t(Delta)
Delta_H <- Delta
# average-linkage
Delta <- as.dist(Delta)
cl <- hclust(Delta, method = "average")
# list of possible clusterings
max.k <- 10  # maximum number of clusters
labs <- matrix(0, nrow = max.k, ncol = n)
for (h in 1:max.k){
  labs[h,] = cutree(cl, k = h)
}
w_avg <- sum(Delta_UT)/(sum(1 - Delta_UT))
losses <- apply(labs, 1, mergeloss, Delta_lt = Delta_H[lower.tri(Delta_H)], w=w_avg)
# default value
min_loss <- which.min(losses)
#c.fold = labs[min_loss,]
# setting equal to number of crescents
c.fold = labs[2,]
# plotting with Binder's Loss
cluster.df <- data.frame(y1 = y[,1],
                         y2 = y[,2],
                         cluster = as.factor(c.VI))
cluster.df %>%
  ggplot(aes(x = y1, y = y2, color = cluster)) + geom_point() + theme_bw() +
  xlab(" ") + ylab(" ") +
  xlim(-2.2,2.2) + ylim(-2.2,2.2) +
  theme(legend.position = "none",
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# elbow plot
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
  ylab("Within Hellinger distance") + xlab("No. of Clusters") +
  scale_x_discrete(limits = as.character(1:max.k)) +
  theme_bw() +
  theme(text = element_text(size=15)) 
# plotting with FOLD
cluster.df$fold <- as.factor(labs[2,])
cluster.df %>%
  ggplot(aes(x = y1, y = y2, color = fold)) + geom_point() + theme_bw() +
  xlab(" ") + ylab(" ") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
# plotting theta
theta_bar = rowMeans(theta, dims = 2)
df.bar <- data.frame(y1 = y[,1], y2=y[,2],
                     t1 = theta_bar[,1], t2 = theta_bar[,2],
                     truth = factor(data$Class),
                     fold = factor(c.fold))
df.bar %>%
  ggplot(aes(x = y1, y = y2)) + geom_point(color = "red") +
  geom_point(aes(x = t1, y = t2), color = "blue")
# plotting localized densities
alpha <- 0.95
rad <- sqrt(qchisq(alpha, df = d))
df.bar$radius <- rep(sqrt(sig_sq)*rad,n) # get 95% density region
colors <- c("Observations" = "black", "Localized Means"="red")
df.bar %>%
  ggplot() + 
  geom_circle(aes(x0 = t1, y0 = t2, r=radius, fill = fold), alpha=0.3,color = "transparent") +
  geom_point(aes(x = t1, y =t2, color = "Localized Means")) +
  geom_point(aes(x = y1, y = y2, color = "Observations"), shape = 1) +
  xlab(" ") + ylab(" ") +
  labs(fill = "Cluster", color = " ") +
  theme_bw() + 
  scale_fill_brewer(palette="Accent") +
  scale_color_manual(values = colors) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        text = element_text(size=15))
# with just clusters
df.bar %>%
  ggplot() + 
  geom_circle(aes(x0 = t1, y0 = t2, r=radius, fill = fold), alpha=0.3,color = "transparent") +
  geom_point(aes(x = t1, y =t2, color = "Localized Means")) +
  xlab(" ") + ylab(" ") +
  theme_bw() + 
  scale_color_manual(values = colors) +
  scale_fill_brewer(palette="Accent") +
  labs(fill = "Cluster", color = " ") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        text = element_text(size=15))
# with just the data
df.bar %>%
  ggplot() +
  geom_point(aes(x = y1, y = y2), shape = 1) + 
  xlab(" ") + ylab(" ") +
  theme_bw() + 
  xlim(-2.2,2.2) + ylim(-2.2,2.2) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=15))
# with just the data and clusters
df.bar %>%
  ggplot() +
  geom_point(aes(x = y1, y = y2, color = truth)) + 
  xlab(" ") + ylab(" ") +
  theme_bw() + 
  xlim(-2.2,2.2) + ylim(-2.2,2.2) + 
  scale_color_brewer(palette="Set1") +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=15))

# over the sth iteration
s <- 1000
theta_s <- theta[,,s]
df.bar$ts1 <- theta_s[,1]
df.bar$ts2 <- theta_s[,2]
df.bar %>%
  ggplot() + 
  geom_circle(aes(x0 = ts1, y0 = ts2, r=radius), alpha=0.01, fill = "grey", color = "transparent") +
  geom_point(aes(x = ts1, y =ts2), color = "red") +
  geom_point(aes(x = y1, y = y2), shape = 1) +
  xlab(" ") + ylab(" ") +
  xlim(-2.2,2.2) + ylim(-2.2,2.2) +
  theme_bw() + 
  scale_fill_brewer(palette="Accent") +
  scale_color_manual(values = colors) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        text = element_text(size=15))
# on average across all samples
df.bar %>%
  ggplot() + 
  geom_circle(aes(x0 = t1, y0 = t2, r=radius), alpha=0.01, fill = "grey", color = "transparent") +
  geom_point(aes(x = t1, y =t2), color = "red") +
  geom_point(aes(x = y1, y = y2), shape = 1) +
  xlab(" ") + ylab(" ") +
  xlim(-2.2,2.2) + ylim(-2.2,2.2) +
  theme_bw() + 
  scale_fill_brewer(palette="Accent") +
  scale_color_manual(values = colors) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        text = element_text(size=15))
# with clusters
df.bar %>%
  ggplot() + 
  geom_circle(aes(x0 = t1, y0 = t2, r=radius, fill = fold), alpha=0.05,color = "transparent") +
  geom_point(aes(x = t1, y =t2, color = "red")) +
  geom_point(aes(x = y1, y = y2), shape = 1) +
  xlab(" ") + ylab(" ") +
  labs(fill = "Cluster", color = " ") +
  theme_bw() + 
  xlim(-2.2,2.2) + ylim(-2.2,2.2) +
  scale_fill_brewer(palette="Accent") +
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        text = element_text(size=15))
# with different allocations
omega_clusters <- matrix(0, nrow = 6, ncol = n)
# hard coding levels of clusters for plot consistency
cfour <- recode(factor(labs[4,]), "4" = "3", "3"="2", "2"="4")
cfive <- recode(factor(labs[5,]), "2" = "4", "5"="3", "4"="2", "3"="5")
csix <- recode(factor(labs[6,]), "5"="2", "6"="3", "4"="5", "3"="4", "2"="6")
omega_clusters <- data.frame(cone = factor(labs[1,]),
                             ctwo = factor(labs[2,]),
                             cthree = factor(labs[3,]),
                             cfour = cfour,
                             cfive = cfive,
                             csix=csix)
customfill <- c("1"="#7FC97F", "2"="#BEAED4",
                "3"="#FDC086", "4"="#FFFF99",
                "5"="#386CB0", "6"="#F0027F",
                "7"="#BF5B17", "8"="#666666") # custom palette
k <- 6
df.bar$cluster <- omega_clusters[,k]
df.bar %>%
  ggplot() + 
  geom_circle(aes(x0 = t1, y0 = t2, r=radius, fill = cluster), alpha=0.05,color = "transparent") +
  geom_point(aes(x = t1, y =t2, color = "red")) +
  geom_point(aes(x = y1, y = y2), shape = 1) +
  xlab(" ") + ylab(" ") +
  labs(fill = "Cluster", color = " ") +
  theme_bw() + 
  xlim(-2.2,2.2) + ylim(-2.2,2.2) +
  scale_fill_manual(values=customfill) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.key.size = unit(1.5, 'cm'),
        text = element_text(size=15))
# credible ball and cGPSM
heights <- rev(cl$height)[1:max.k]
w_vals <- heights/(1-heights)
cgsamps <- getcGsamps_loc(theta=theta, w=w_vals[3], d=2, max.k=10, Sigma=Sigma) # to match figure in main article
# cgsamps <- read.csv("output/cells/cgsamps.csv")
cb <- cGball(labs[3,],cgsamps)
foldpsm <- cGpsm(cgsamps)
# plot PSM
c.fold.ordered <- order(labs[3,])
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
  y,
  cfold = factor(labs[3,]),
  choriz = factor(cb$c.horiz[1,]),
  c.upperv = factor(cb$c.uppervert[1,]),
  c.lowerv = factor(cb$c.lowervert[1,]))
foldcl <- ggplot(ball.df, aes(x = X1, y = X2, color = cfold)) + geom_point() + 
  labs(title=" ", color = "Cluster", caption = latex2exp::TeX("$\\textbf{c}_{\\textrm{FOLD}}$.")) +
  xlab(" ") + ylab(" ") +
  theme_bw() + 
  #scale_color_manual(values=customfill) +
  scale_color_brewer(palette="Dark2")+
  theme(axis.line = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        text = element_text(size=14),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.caption = element_text(hjust=0.5))
horizcl <- ggplot(ball.df, aes(x = X1, y = X2, color = choriz)) + geom_point() + 
  labs(caption="95% credible ball horizontal bound.", color = "Cluster", title = " ") + 
xlab(" ") + ylab(" ") +
  scale_color_brewer(palette="Dark2")+
  theme_bw() + 
  theme(axis.line = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        text = element_text(size=14),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.caption = element_text(hjust=0.5))
uppervcl <- ggplot(ball.df, aes(x = X1, y = X2, color = c.upperv)) + geom_point() + 
  labs(caption="95% credible ball vertical upper bound.", color = "Cluster", title = " ") + 
xlab(" ") + ylab(" ") +
  scale_color_brewer(palette="Dark2")+
  theme_bw() + 
  theme(axis.line = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        text = element_text(size=14),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.caption = element_text(hjust=0.5))
lowervcl <- ggplot(ball.df, aes(x = X1, y = X2, color = c.lowerv)) + geom_point() + 
  labs(caption="95% credible ball vertical lower bound.", color = "Cluster", title = " ") +
xlab(" ") + ylab(" ") +
  theme_bw() + 
  scale_color_brewer(palette="Dark2")+
  theme(axis.line = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        text = element_text(size=14),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.caption = element_text(hjust=0.5))
# plotting point estimate with credible ball 
plot_grid(foldcl, horizcl, lowervcl, uppervcl, labels = c("(a)", "(b)", "(c)", "(d)"))
