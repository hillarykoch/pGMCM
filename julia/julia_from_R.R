library(JuliaCall)
library(mvtnorm)
source("~/Box Sync/School/research - Qunhua/Project_2/pGMCM/R/R/ptgibbs.R")

# Prepare some inputs
dat <- data.frame(
    rbind(rmvnorm(100, mean = c(0,0), sigma = diag(2)),
          rmvnorm(200, mean = c(4,4), sigma = diag(2) + .7))
)

kappa0 <- nu0 <- round(c(1/3,2/3)*nrow(dat))
mu0 <- rbind(c(0,0), c(4,4))
Psi0 <- base::array(data = c(2, 0, 0, 2, 1.7, .7, .7, 1.7), dim = c(2, 2, 2))
alpha <- c(1/3,2/3)
nw <- 2
nt <- 3
nstep <- 50
burnin <- 10


# Test out setting up julia and running the chain
prepare_julia()
chain <- run_ptgibbs(dat, kappa0, mu0, Psi0, alpha, nw, nt, nstep, burnin)

# Process the results
get_mu_chain(chain, 1, 1)
