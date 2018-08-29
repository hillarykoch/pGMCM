# Do all pairwise fits on 100 replicates of the simulation
library(magrittr)
library(abind)
library(pGMCM)
library(RColorBrewer)
library(foreach)
library(parallel)
library(doParallel)
set.seed(719)
pal <- c(brewer.pal(8, "Pastel1"), brewer.pal(8, "Dark2"))

# Choose which subset of indices are the truth. Here, no negative associations
n <- 3000
d <- 5
c_vis <- expand.grid(rep(list(-1:1), d)) %>% as.matrix
l <- nrow(c_vis)

neg <- apply(c_vis, 1, function(X)
    all(X != -1))
nonnegidx <- which(neg)

# trueidx are the indices from c_vis which correspond to the true underlying
# associations
trueidx <- sample(nonnegidx, size = 15, replace = FALSE)
prop <- rep(0, l) %>% replace(trueidx, 1)

mu <- 2.5
sigma <- 0.8
rho <- 0.55

replicates <- 100
sim <- replicate(replicates,
                 rconstr_GMCM(n, prop, mu, sigma, rho, d),
                 simplify = FALSE)

# Quick visualization
par(mfrow = c(2, 2))
plot(sim$z[, c(1, 2)], col = pal[sim$cluster], pch = 16)
plot(sim$z[, c(1, 3)], col = pal[sim$cluster], pch = 16)
plot(sim$z[, c(1, 4)], col = pal[sim$cluster], pch = 16)
plot(sim$z[, c(1, 5)], col = pal[sim$cluster], pch = 16)

# setup parallel backend to use given number of cores
cores <- howmany
cluster <- makeCluster(cores)
registerDoParallel(cluster)
clusterEvalQ(cluster, library(pGMCM))
pairs <- combn(d, 2)

fits <-
    foreach(j = seq(replicates), .packages = c("foreach")) %dopar% {
        foreach(l = seq(ncol(pairs))) %do% {
            fconstr_pGMCM(
                sim[[j]]$z[, pairs[, l]],
                tol = 1e-04,
                stepmax = 50,
                itermax = 100,
                convCrit = "GMM",
                penaltyType = "SCAD"
            )
        }
    }

stopCluster(cluster)

save(fits, file = "naive_fits.Rdata")
