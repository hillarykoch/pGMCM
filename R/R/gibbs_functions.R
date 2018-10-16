# Functions needed to process process and run gibbs sampler

# Get expected value of prior probability on full mixing proportions
# Based on empirical concordance with various paths in red_class
get_prior_prop <- function(red_class, fits, d, n, dist_tol = 0){
    mus <- map(fits, "mu")
    labels <- map(fits, "cluster") %>% simplify2array

    prior_count <- cget_prior_count(red_class, mus, labels, d, n, dist_tol)
    (prior_count/n)/(sum(prior_count/n))
}

# Get expected value of prior probability on full mixing proportions
get_prior_prop_OLD <- function(fits, red_class, d) {
    mus <- purrr::map(fits, "mu")
    props <- purrr::map(fits, "prop")
    cget_prior_prop(red_class, mus, props, d)
}

# Averaging mus, sigmas in each dimension, recording rhos directly
# The covariance is NOT handled correctly at the moment
# Maybe put this is c++ later
get_hyperparams <- function(fits, d, red_class) {
    mu0_temp <- rep(NA, d)
    Psi0_temp <- matrix(rep(NA, d^2), nrow = d, ncol = d)
    spl <- names(fits) %>%
        str_split(pattern = "_") %>%
        purrr::map(as.numeric) %>%
        abind(along = 2) %>%
        t
    mus <- purrr::map(fits, "mu") %>% purrr::map(abs)
    sigmas <- purrr::map(fits, "sigma")
    rhos <- purrr::map(fits, "rho") %>% purrr::map(abs)

    # For each dimension
    for(i in seq(d)) {
        idx <- t(apply(spl, 1, function(X) X == i))
        muvec <- sigmavec <- rep(NA, sum(idx))

        # In each pairwise fit that models that dimension
        count <- 1
        for(j in 1:nrow(idx)) {
            if(any(idx[j,])) {
                submu <- mus[[j]][,which(idx[j,])]
                subsigma <- sigmas[[j]][,which(idx[j,])]

                muvec[count] <- submu[submu != 0][1]
                sigmavec[count] <- subsigma[subsigma != 1][1]
                count <- count+1
            }
        }
        mu0_temp[i] <- mean(muvec)
        Psi0_temp[i,i] <- mean(sigmavec) # This covariance probably isnt correct
    }

    for(i in seq(nrow(spl))) {
        # This covariance probably isnt correct
        Psi0_temp[spl[i,1], spl[i,2]] <- Psi0_temp[spl[i,2], spl[i,1]] <- rhos[[i]]
    }


    mu0 <- matrix(rep(0, length(red_class)),
                  nrow = nrow(red_class),
                  ncol = ncol(red_class))
    Psi0 <- base::array(rep(0, length(Psi0_temp)*nrow(red_class)),
                        dim = c(d, d, nrow(red_class)))

    # Update mu0 based on association patterns
    for(i in seq(ncol(red_class))) {
        mu0[red_class[,i] == 1, i] <- mu0_temp[i]
        mu0[red_class[,i] == -1, i] <- -mu0_temp[i]
    }

    # Update Psi0 based on association patterns
    for(i in seq(nrow(red_class))) {
        diag(Psi0[,,i]) <- diag(Psi0_temp)*red_class[i,] %>%
            abs
        diag(Psi0[,,i]) <- replace(diag(Psi0[,,i]),
                                   list = diag(Psi0[,,i]) == 0,
                                   values = 1)

        # Don't update rho if there arent multiple associations
        assoc_idx <- red_class[i,] != 0
        if(sum(assoc_idx) > 1) {
            wai <- which(assoc_idx)

            # for every pairwise associated combination
            cmb <- combn(wai,2)
            for(j in seq(ncol(cmb))) {
                Psi0[cmb[1, j], cmb[2, j], i] <-
                    Psi0[cmb[2, j], cmb[1, j], i] <-
                    red_class[i,cmb[1,j]] * red_class[i,cmb[2,j]] * Psi0_temp[cmb[1, j], cmb[2, j]]
            }
        }
    }
    list("mu0" = mu0, "Psi0" = Psi0)
}

get_kappa <- function(z, nclass) {
    kappa <- rep(0, nclass)
    runlens <- rle(sort(z))
    kappa[runlens$values] <- runlens$lengths
    kappa
}

draw_mix_prop <- function(alpha, z){
    # Counting the number of observations that belong to each class
    runlens <- rle(sort(z))
    zcoll <- rep(0, length(alpha))
    zcoll[runlens$values] <- runlens$lengths

    # Sampling from dirichlet dist with param alpha + zcoll
    az <- alpha + zcoll
    nimble::rdirch(n=1, alpha = az/sum(az))
}

# Mean and variance
draw_NIW <- function(x, hyp, z) {
    nz <- hyp$kappa#get_kappa(z, nclass)
    xbar <- sapply(seq(ncol(x)), function(X) tapply(x[,X], z, mean))

    lapply(1:nrow(hyp$hyp$mu0), function(X)
        LaplacesDemon::rnorminvwishart(mu0 = (hyp$kappa[X]*hyp$hyp$mu0[X,] + nz[X]*xbar[X,])/(hyp$kappa[X] + nz[X]),
                                       lambda = hyp$kappa[X] + nz[X],
                                       nu = hyp$kappa[X] + nz[X],
                                       S = hyp$hyp$Psi0[,,X]*hyp$kappa[X]))
}

updatez <- function(data, NIW) {
    d <- sapply(seq_along(NIW),
                function(X) mvtnorm::dmvnorm(data,
                                             mean = NIW[[X]]$mu,
                                             sigma = NIW[[X]]$Sigma))
    apply(d, 1, which.max)
}

# run the gibbs sampler
run_gibbs <- function(x, prior_prop, fits, d, red_class, nsamp) {
    n <- nrow(x)
    nclass <- nrow(red_class)
    z_init <- sample(seq(nclass), replace = TRUE, size = n, prob = prior_prop)
    kappa <- get_kappa(z_init, nclass) # kappa = nu. get_kappa can also be used to update n^z
    hyp <- list("kappa" = kappa, "hyp" = get_hyperparams(fits, d, red_class))
    prop0 <- draw_mix_prop(alpha = prior_prop, z = z_init) # can probably update this to just pass zcoll
    NIW <- draw_NIW(x, hyp, z_init)
    z <- updatez(x, NIW)

    # Run the gibbs sampler (with a progress bar!)
    param.tr <- list()
    param.tr[[1]] <- list("z" = z, "mix_prop" = prop0, "NIW" = NIW)

    pb <- txtProgressBar(min = 2, max = nsamp, style = 3)
    for(i in 2:nsamp){
        mix_prop <- draw_mix_prop(prior_prop, param.tr[[i-1]]$z)
        NIW <- draw_NIW(x, hyp, param.tr[[i-1]]$z)
        zstar <- updatez(x, NIW)
        param.tr[[i]] <- list("z" = zstar, "mix_prop" = mix_prop, "NIW" = NIW)
        setTxtProgressBar(pb,i)
    }
    close(pb)

    param.tr
}

