library(pGMCM)
library(abind)
library(tidyverse)

## generate data for 3-component mixture
set.seed(123)
n <- 1500    # sample size
m <- 3      # number of components
d <- 2      # dimension

prop <- c(1/3, 1/3, 1/3)
mu <- list(c(-1, 1), c(1, 1), c(0, -sqrt(2)))
sigma1 <- rbind(c(0.65, 0.7794), c(0.7794, 1.55))
sigma2 <- rbind(c(0.65, -0.7794), c(-0.7794, 1.55))
sigma3 <- diag(c(2, 0.2))
sigma <- list(sigma1, sigma2, sigma3)

sim <- rGMCM(n, prop, mu, sigma)
data <- sim$data
clust <- sim$component
u <- sim$u   # rescaled empirical marginal distribution
z <- sim$z   # pseudo data

x <- as.matrix(sim$data)
k <- 5
out <- fpGMCM(x, k)
out2 <- fpGMM(x, k)
#out3 <- fit.full.GMCM(x,3, method = "PEM")

pdf(file = "../triangle.pdf")
par(mfrow = c(2,2))
pal <- get_pals(6)
plot(x, col = pal[sim$cluster+1], xlab = expression(x[1]), ylab = expression(x[2]), main = "Observed data", pch = 16, cex = .7)
plot(u, col = pal[sim$cluster+1], xlab = expression(x[1]), ylab = expression(x[2]), main = "Copula process", pch = 16, cex = .7)
plot(x, col = pal[(out2$cluster+1)], xlab = expression(x[1]), ylab = expression(x[2]), main = "pGMM classification", pch = 16, cex = .7)
plot(x, col = pal[out$cluster+1], xlab = expression(x[1]), ylab = expression(x[2]), main = "pGMCM classification", pch = 16, cex = .7)
dev.off()

# generate data for 9 component mixture
set.seed(123)
n <- 2700
k <- 9
d <- 2

prop <- rep(1/9,9)
mu <- 3
sigma <- 1.1
rho <- .7

sim <- rconstr_GMCM(n, prop, mu, sigma, rho, d)
data <- sim$data
clust <- sim$component
u <- sim$u
z <- sim$z

x <- as.matrix(sim$data)
k <- 9
lambda <- seq(.281,5,length.out = 5)

out <- fconstr_pGMCM(x, k, lambda, tol = 1e-04, stepmax = 10, itermax = 50)
out2 <- fconstr_pGMM(x, k, lambda, tol = 1e-04, itermax = 100)

pdf(file = "../nine.pdf")
pal <- get_pals(5)
par(mfrow = c(2,2))
plot(x, col = pal[sim$cluster+1], xlab = expression(x[1]), ylab = expression(x[2]), main = "Observed data", pch = 16, cex = .7)
plot(u, col = pal[sim$cluster+1], xlab = expression(x[1]), ylab = expression(x[2]), main = "Copula process", pch = 16, cex = .7)
plot(x, col = pal[(out2$cluster+1)], xlab = expression(x[1]), ylab = expression(x[2]), main = "pGMM classification", pch = 16, cex = .7)
plot(x, col = pal[out$cluster+1], xlab = expression(x[1]), ylab = expression(x[2]), main = "pGMCM classification", pch = 16, cex = .7)
dev.off()


# remove some data from 9 component mixture
set.seed(123)
n <- 1200
k <- 9
d <- 2

prop <- c(1,0,1,0,1,0,0,1,0)
mu <- 3
sigma <- 1.1
rho <- .7

sim <- rconstr_GMCM(n, prop, mu, sigma, rho, d)
data <- sim$data
clust <- sim$component
u <- sim$u
z <- sim$z

x <- as.matrix(sim$data)
lambda <- seq(.281,5,length.out = 5)

out <- fconstr_pGMCM(x, lambda, tol = 1e-04, stepmax = 10, itermax = 50)
out2 <- fconstr_pGMM(x, lambda, tol = 1e-04, itermax = 100)

pdf(file = "../subnine.pdf")
pal <- get_pals(1)
par(mfrow = c(2,2))
plot(x, col = pal[sim$cluster+2], xlab = expression(x[1]), ylab = expression(x[2]), main = "Observed data", pch = 16, cex = .7)
plot(u, col = pal[sim$cluster+2], xlab = expression(x[1]), ylab = expression(x[2]), main = "Copula process", pch = 16, cex = .7)
plot(x, col = pal[(out2$cluster+2)], xlab = expression(x[1]), ylab = expression(x[2]), main = "pGMM classification", pch = 16, cex = .7)
plot(x, col = pal[out$cluster+2], xlab = expression(x[1]), ylab = expression(x[2]), main = "pGMCM classification", pch = 16, cex = .7)
dev.off()
