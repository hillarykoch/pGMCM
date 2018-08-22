library(pGMCM)
library(magrittr)
library(abind)
set.seed(637)
pal <- sample(get_pals(5), 9, replace = F)

sim <- rconstr_GMCM(2000, prop = rep(1/9,9), mu = 2, sigma = .9, rho = .75, d=2)
x <- sim$data

# Fit with both stopping rules
fit1 <- fconstr_pGMCM(x, lambda = 0, tol = 1e-03, stepmax = 20, itermax = 70, convCrit = "GMCM")
fit2 <- fconstr_pGMCM(x, lambda = 0, tol = 1e-03, stepmax = 20, itermax = 70, convCrit = "GMM")

par(mfrow = c(1,3))
plot(x, col = pal[sim$cluster], main = "true clustering")
plot(x, col = pal[fit1$cluster], main = "GMCM stopping rule clustering")
plot(x, col = pal[fit2$cluster], main = "GMM stopping rule clustering")

# GMCM stopping rule performs slightly better
sum((fit1$cluster-sim$cluster)==0)
sum((fit2$cluster-sim$cluster)==0)

# Test fconstr0 -- can't tell if this is my bug or just an unidentifiability issue
fit3 <- fconstr0_pGMCM(x, lambda = 0, tol = 1e-03, stepmax = 25, itermax = 50, convCrit = "GMCM")
fit4 <- fconstr0_pGMCM(x, lambda = 0, tol = 1e-03, stepmax = 25, itermax = 50, convCrit = "GMM")

plot(x, col = pal[sim$cluster], main = "true clustering")
plot(x, col = pal[fit3$cluster], main = "GMCM stopping rule clustering")
plot(x, col = pal[fit4$cluster], main = "GMM stopping rule clustering")
