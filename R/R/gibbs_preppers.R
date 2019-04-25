# Functions needed to process process and run gibbs sampler

# Get expected value of prior probability on full mixing proportions
# Based on empirical concordance with various paths in red_class
get_prior_prop <- function(red_class, fits, d, n, dist_tol = 0, MAP = TRUE) {
    mus <- map(fits, "mu")
    
    if (MAP) {
        labels <- purrr::map(fits, "cluster") %>% simplify2array    
    } else {
        labels <- purrr::map(fits, "post_prob") %>%
            lapply(function(X) apply(X, 1, function(Y) base::sample(1:ncol(X), size = 1, prob = Y))) %>%
            simplify2array
    }
    
    prior_count <-
        cget_prior_count(red_class, mus, labels, d, n, dist_tol)
    (prior_count / n) / (sum(prior_count / n))
}

# This function will filter out classes that are really similar to each other
# since they probably are not that informative and this task of computing
# prior proportions (alpha) is super burdensome when there are too many classes
reduce_by_hamming <- function(red_class, hamming_tol = 1, force_canonical = TRUE) {
    M <- nrow(red_class)
    fullidx <- seq(M)
    out <- rep(NA, M)
    
    count <- 1
    while(length(fullidx) != 0) {
        candidate_indices <- creduce_by_hamming(as.matrix(red_class), fullidx-1, hamming_tol, M)
        candidate_indices <- fullidx[which(candidate_indices == 1)]
        
        if(length(candidate_indices) == 1) {
            keepidx <- candidate_indices
        } else {
            keepidx <- sample(candidate_indices, size = 1)    
        }
        
        out[count] <- keepidx
        fullidx <- fullidx[!(fullidx %in% candidate_indices)]
        count <- count + 1
    }
    
    outidx <- out[!is.na(out)]
    
    # If we have canonical behavior, keep it
    if(force_canonical) {
        allneg <- which(sapply(1:M, function(X) all(red_class[X,] == -1)))
        allone <- which(sapply(1:M, function(X) all(red_class[X,] == 0)))
        allpos <- which(sapply(1:M, function(X) all(red_class[X,] == 1)))
        outidx <- unique(c(outidx, allneg, allpos, allone))
    }
    
    red_class[outidx,]
}


# Averaging mus, sigmas in each dimension, recording rhos directly
# The covariance is NOT handled correctly at the moment
# Maybe put this is c++ later
get_hyperparams <- function(fits, d, red_class) {
    mu0_temp <- rep(NA, d)
    Psi0_temp <- matrix(rep(NA, d ^ 2), nrow = d, ncol = d)
    spl <- names(fits) %>%
        str_split(pattern = "_") %>%
        purrr::map(as.numeric) %>%
        abind::abind(along = 2) %>%
        t
    mus <- purrr::map(fits, "mu") %>% purrr::map(abs)
    sigmas <- purrr::map(fits, "sigma")
    rhos <- purrr::map(fits, "rho") %>% purrr::map(abs)
    
    # For each dimension
    for (i in seq(d)) {
        idx <- t(apply(spl, 1, function(X)
            X == i))
        muvec <- sigmavec <- rep(NA, sum(idx))
        
        # In each pairwise fit that models that dimension
        count <- 1
        for (j in 1:nrow(idx)) {
            if (any(idx[j, ])) {
                submu <- mus[[j]][, which(idx[j, ])]
                subsigma <- sigmas[[j]][, which(idx[j, ])]
                
                muvec[count] <- submu[submu != 0][1]
                sigmavec[count] <- subsigma[subsigma != 1][1]
                count <- count + 1
            }
        }
        mu0_temp[i] <- mean(muvec)
        Psi0_temp[i, i] <-
            #mean(sigmavec) # This covariance probably isnt correct
            max(sigmavec) * 5 # This covariance probably isnt correct
    }
    
    for (i in seq(nrow(spl))) {
        # This covariance probably isnt correct
        Psi0_temp[spl[i, 1], spl[i, 2]] <-
            Psi0_temp[spl[i, 2], spl[i, 1]] <- rhos[[i]] #/ 10
    }
    
    
    mu0 <- matrix(rep(0, length(red_class)),
                  nrow = nrow(red_class),
                  ncol = ncol(red_class))
    Psi0 <- base::array(rep(0, length(Psi0_temp) * nrow(red_class)),
                        dim = c(d, d, nrow(red_class)))
    
    # Update mu0 based on association patterns
    for (i in seq(ncol(red_class))) {
        mu0[red_class[, i] == 1, i] <- mu0_temp[i]
        mu0[red_class[, i] == -1, i] <- -mu0_temp[i]
    }
    
    # Update Psi0 based on association patterns
    for (i in seq(nrow(red_class))) {
        diag(Psi0[, , i]) <- as.numeric(unname(diag(Psi0_temp) * red_class[i, ] %>%
            abs))
        diag(Psi0[, , i]) <- replace(diag(Psi0[, , i]),
                                     list = diag(Psi0[, , i]) == 0,
                                     values = 1)
        
        # Don't update rho if there arent multiple associations
        assoc_idx <- red_class[i, ] != 0
        if (sum(assoc_idx) > 1) {
            wai <- which(assoc_idx)
            
            # for every pairwise associated combination
            cmb <- combn(wai, 2)
            for (j in seq(ncol(cmb))) {
                Psi0[cmb[1, j], cmb[2, j], i] <-
                    Psi0[cmb[2, j], cmb[1, j], i] <-
                    red_class[i, cmb[1, j]] * red_class[i, cmb[2, j]] * Psi0_temp[cmb[1, j], cmb[2, j]]
            }
        }
    }
    list("mu0" = mu0, "Psi0" = Psi0)
}

# Count how many observations belong to each class
get_kappa <- function(z, nclass) {
    kappa <- rep(0, nclass)
    runlens <- rle(sort(z))
    kappa[runlens$values] <- runlens$lengths
    kappa
}