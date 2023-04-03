require(mvtnorm)
require(sn)
require(ggplot2)
require(cowplot)
# normals
d.gaussmix <- function(x, mu, Sigma, Pi){
  fx <- Pi[1] * mvtnorm::dmvnorm(x=x, mean = mu[1,], sigma = Sigma[,,1]) +
    Pi[1] * mvtnorm::dmvnorm(x=x, mean = mu[2,], sigma = Sigma[,,2]) +
    Pi[1] * mvtnorm::dmvnorm(x=x, mean = mu[3,], sigma = Sigma[,,3]) 
  return(fx)
}
# parameters
mu <- matrix(c(6.5,5,
               0,0,
               -5,-5), byrow = T, nrow = 3)
Sigma <- array(c(1,0,0,1,5,0,0,2,3,0,0,1), c(2,2,3))
Pi <- c(0.45, 0.25, 0.3)
# plot
x1 <- seq(-10,10, by = 0.05)
x2 <- seq(-10,10, by = 0.05)
xgrid <- expand.grid(x1,x2)
x3 <- apply(xgrid, 1, d.gaussmix, 
            mu=mu,
            Sigma=Sigma,
            Pi=Pi)
xgrid <- cbind(xgrid, x3)
colnames(xgrid) <- c("x1", "x2", "x3")
normal_contours = ggplot(xgrid, aes(x = x1, y = x2, z = x3)) + geom_contour() +
  labs(title = "Contours of Skew-Symmetric Mixture") +
  xlab(" ") + ylab(" ") +
  theme(text = element_text(size = 14),
        title = element_text(size = 16)) +
  theme_bw()


# skew normals
x1 <- seq(-10,10, by = 0.01)
x2 <- seq(-9,9, by = 0.01)
dmsn_mix <- function(x, Pi, xi, Omega, alpha) {
  x3 <- Pi[1]* dmsn(x = x, xi = xi[1,], Omega = Omega[,,1], alpha = alpha[1,]) +
    Pi[2]* dmsn(x = x, xi = xi[2,], Omega = Omega[,,2], alpha = alpha[2,]) +
    Pi[3] * dmsn(x = x, xi = xi[3,], Omega = Omega[,,3], alpha = alpha[3,])
  return(x3)
}
xgrid <- expand.grid(x1,x2)
x3 <- apply(xgrid, 1, dmsn_mix, 
            Pi = c(0.45, 0.25, 0.3),
            xi = matrix(c(6.5,5,
                          0,0,
                          -5,-5), byrow = T, nrow = 3),
            Omega = array(c(1,0,0,1,5,0,0,2,3,0,0,1), c(2,2,3)),
            alpha = matrix(c(1,1,
                             -10,15,
                             4,-17), byrow = T, nrow = 3))
xgrid <- cbind(xgrid, x3)
colnames(xgrid) <- c("x1", "x2", "x3")
skewnormal_contours = ggplot(xgrid, aes(x = x1, y = x2, z = x3)) + geom_contour() +
  labs(title = "Contours of Skew Normal Mixture") + 
  xlab(" ") + ylab(" ") +
  theme(text = element_text(size = 14),
        title = element_text(size = 16)) +
  theme_bw()

# skew-symmetric
d.multivar.gmm <- function(x, p, Mu, Sigma){
  s <- sample(1:2, size = n, replace = T, prob = c(p,1-p))
  x <- p * mvtnorm::dmvnorm(x = x, mean = Mu[1,], sigma = Sigma[,,1]) +
    (1-p) * mvtnorm::dmvnorm(x=x, mean = Mu[2,], sigma = Sigma[,,2])
  return(x)
}
d.nongauss <- function(x, Pi, xi, Omega, alpha, Mu_1, Sigma_1, Mu, Sigma, p){
  # omega = vector of probabilties for mixture
  # p_mix, mu_mix, sigma_mix = parameters for 3 component mixture
  # mu, sigma = normal parameters
  # shape, rate = gamma parameters
  require(sn)
  x <- Pi[1] * dmsn(x=x, xi = xi[1,], Omega = Omega[,,1], alpha = alpha[1,]) +
    Pi[2] * dmsn(x=x, xi = xi[2,], Omega = Omega[,,2], alpha = alpha[2,]) +
    Pi[3] * mvtnorm::dmvnorm(x=x, mean = Mu_1, sigma = Sigma_1) +
    Pi[4] * d.multivar.gmm(x = x, p = p, Mu = Mu, Sigma = Sigma)
  return(x=x)
}
x1 <- seq(-7,7, by = 0.01)
x2 <- seq(-7,7, by = 0.01)
xgrid <- expand.grid(x1,x2)
x3 <- apply(xgrid, 1, d.nongauss, 
            Pi = c(0.2, 0.3, 0.15, 0.35),
            xi = matrix(c(2.5,3.5,
                          0,-3.5), byrow = T, nrow = 2),
            Omega = array(c(1,0,0,1,5,0,0,2), c(2,2,2)),
            alpha = matrix(c(-10,15,
                             4,-17), byrow = T, nrow = 2),
            Mu_1 = c(-4,2.5),
            Sigma_1 = matrix(c(1/2, 1/2, 1/2, 5/2), byrow = T, nrow = 2),
            Mu = matrix(c(2.325, 1.085, 4.381, 2.009), nrow = 2),
            Sigma = array(c(0.2, 0, 0, 0.8, 0.7, 0, 0, 0.6), c(2,2,2)),
            p = 1/3)
xgrid <- cbind(xgrid, x3)
colnames(xgrid) <- c("x1", "x2", "x3")
skewsymm_contours = ggplot(xgrid, aes(x = x1, y = x2, z = x3)) + geom_contour() +
  labs(title = "Contours of Skew-Symmetric Mixture") +
  xlab(" ") + ylab(" ") +
  theme(text = element_text(size = 14),
        title = element_text(size = 16)) +
  theme_bw()
plot_grid(normal_contours, skewnormal_contours, skewsymm_contours, ncol = 1, labels = c("(a)", "(b)", "(c)"))