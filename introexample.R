# packages
library(foldcluster)
library(ggplot2)
library(dplyr)
library(mcclust.ext)
library(RColorBrewer)
library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)
# functions
# sample Multivariate Normal
multivar.gmm <- function(n, p, Mu, Sigma){
  s <- sample(1:2, size = n, replace = T, prob = c(p,1-p))
  x <- (s==1) * rmvnorm(n = n, mean = Mu[1,], sigma = Sigma[,,1]) +
    (s==2) * rmvnorm(n=n, mean = Mu[2,], sigma = Sigma[,,2])
  return(x)
}
nongauss <- function(n, Pi, xi, Omega, alpha, Mu_1, Sigma_1, Mu, Sigma, p){
  # omega = vector of probabilties for mixture
  # p_mix, mu_mix, sigma_mix = parameters for 3 component mixture
  # mu, sigma = normal parameters
  # shape, rate = gamma parameters
  require(sn)
  K <- 4
  s <- sample(1:K, size = n, replace = T, prob = Pi)
  x <- (s==1) * rmsn(n=n, xi = xi[1,], Omega = Omega[,,1], alpha = alpha[1,]) +
    (s==2) * rmsn(n=n, xi = xi[2,], Omega = Omega[,,2], alpha = alpha[2,]) +
    (s==3) * rmvnorm(n=n, mean = Mu_1, sigma = Sigma_1) +
    (s==4) * multivar.gmm(n = n, p = p, Mu = Mu, Sigma = Sigma)
  return(list(data = x, s = s))
}
# sampling
set.seed(1996)
samp <- nongauss(n = 2000,
                 Pi = c(0.2, 0.3, 0.15, 0.35),
                 xi = matrix(c(2.5,3.5,
                               0,-3.5), byrow = T, nrow = 2),
                 Omega = array(c(1,0,0,1,5,0,0,2), c(2,2,2)),
                 alpha = matrix(c(-10,15,
                                  4,-17), byrow = T, nrow = 2),
                 Mu_1 = c(-4,2.5),
                 Sigma_1 = matrix(c(1/2, 1/2, 1/2, 5/2), byrow = T, nrow = 2),
                 Mu = matrix(rnorm(4, mean = 1, sd = 1), nrow = 2),
                 Sigma = array(c(0.2, 0, 0, 0.8, 0.7, 0, 0, 0.6), c(2,2,2)),
                 p = 1/3)
s <- samp$s
y <- samp$data
y <- scale(y)
n <- nrow(y)
d <- 2
# fitting a GMM
S <- 25000 # iterations
B <- 1000 # burnin
L <- 10
# fitting
fit <- mvnorm_gibbs(S = S,
                    y = y,
                    L = L,
                    alpha = rep(1/L,L), # rep((d/2)+2,L)
                    w = rep(0,d),
                    kappa = 1,
                    r = d+2,
                    C = diag(1,d),
                    stops = 50)
# burn-out
z <- fit$z[-(1:B),]
# thin
trip_ind <- seq(4,nrow(z),by=4)
z <- z[trip_ind,]
# cluster
# W & G estimate -- VI
c.psm <- mcclust::comp.psm(z)
c.mv <- mcclust.ext::minVI(psm = c.psm,
                           cls.draw = z,
                           max.k = 10)
c.VI <- c.mv$cl
# W & G estimate -- Binder's
c.mv <- mcclust::minbinder(psm = c.psm,
                           cls.draw = z,
                           max.k = 10)
c.B <- c.mv$cl
# plotting with Binder's Loss
y_raw <- samp$data
cluster.df <- data.frame(y1 = y_raw[,1],
                         y2 = y_raw[,2],
                         cluster = as.factor(c.B))
cluster.df %>%
  ggplot(aes(x = y1, y = y2, color = cluster)) + geom_point() + theme_bw() +
  xlab(" ") + ylab(" ") +
  theme(legend.position = "none",
        #axis.text = element_blank(),
        #axis.ticks = element_blank(),
        #panel.border = element_blank(),
        text = element_text(size=13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
