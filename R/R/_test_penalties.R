library(pGMCM)
library(magrittr)
library(abind)
set.seed(637)
pal <- sample(get_pals(5), 9, replace = F)

sim <- rconstr_GMCM(2000, prop = c(0,0,0,0,1,1,0,1,1), mu = 2, sigma = .9, rho = .75, d=2)
x <- sim$z

# Fit with both penalty types
fit1 <- fconstr_pGMCM(x, lambda = NULL, tol = 1e-04,
                      stepmax = 20, itermax = 100,
                      convCrit = "GMCM", penaltyType = "SCAD",
                      trace_params = TRUE)
fit2 <- fconstr_pGMCM(x, lambda = NULL, tol = 1e-04,
                      stepmax = 20, itermax = 100,
                      convCrit = "GMCM", penaltyType = "LASSO",
                      trace_params = TRUE)

par(mfrow = c(1,3))
plot(x, col = pal[sim$cluster], main = "true clustering")
plot(x, col = pal[fit1$cluster], main = "SCAD penalty clustering")
plot(x, col = pal[fit2$cluster], main = "LASSO penalty clustering")

# penalties perform equally in misclassification rates
sum((fit1$cluster-sim$cluster)==0)
sum((fit2$cluster-sim$cluster)==0)
