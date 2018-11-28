propose_prop <- function(prop,
                         tune_var = 2) {
    log_prop <- log(prop)
    log_proposal <- rmvnorm(1, mean = log_prop, sigma = diag(sqrt(tune_var), length(log_prop)))
        # sapply(prop, function(X)
        #     rnorm(1, mean = X, sd = tune_var))
    proposal <- exp(log_proposal)
    proposal / sum(proposal)
}

# Automatic tuning variance
update_var <- function(cur_var, acpt_rt, opt_rt = .3, gamma1) {
    exp(log(cur_var) + gamma1 * (acpt_rt - opt_rt))
}

# Gibbs update for NIW
draw_NIW2 <- function(x, hyp, z) {
    kappa <- hyp$kappa
    nclass <- length(kappa)
    nz <- get_kappa(z, nclass)
    xbar <- sapply(seq(ncol(x)), function(X)
        tapply(x[, X], z, mean))
    
    # If there are enough (D) observations in the class, sample from the posterior
    # Otherwise, just draw from the prior
    lapply(seq(nrow(hyp$hyp$mu0)), function(X)
        # If cluster X has enough members to evaluate
        if (nz[X] >= ncol(x)) {
            # If cluster X doesn't hold all the observations
            if(!is.null(dim(xbar))) {
                LaplacesDemon::rnorminvwishart(
                    mu0 = (kappa[X] * hyp$hyp$mu0[X, ] + nz[X] * xbar[as.character(X), ]) / (kappa[X] + nz[X]),
                    lambda = hyp$kappa[X] + nz[X],
                    nu = hyp$kappa[X] + nz[X],
                    S = hyp$hyp$Psi0[, , X] * kappa[X] +
                        t(as.matrix(x)[z==X,] - xbar[as.character(X),]) %*% (as.matrix(x)[z==X,] - xbar[as.character(X),]) +
                        (kappa[X] * nz[X]) / (kappa[X] + nz[X]) *
                        outer(xbar[as.character(X),] - hyp$hyp$mu0[X,],
                              xbar[as.character(X),] - hyp$hyp$mu0[X,])
                )
            } else {
                LaplacesDemon::rnorminvwishart(
                    mu0 = (kappa[X] * hyp$hyp$mu0[X, ] + nz[X] * xbar) / (kappa[X] + nz[X]),
                    lambda = hyp$kappa[X] + nz[X],
                    nu = hyp$kappa[X] + nz[X],
                    S = hyp$hyp$Psi0[, , X] * kappa[X] +
                        t(as.matrix(x)[z==X,] - xbar) %*% (as.matrix(x)[z==X,] - xbar) +
                        (kappa[X] * nz[X]) / (kappa[X] + nz[X]) *
                        outer(xbar - hyp$hyp$mu0[X,], xbar - hyp$hyp$mu0[X,])
                )
            }
        } else {
            # Drawwing from the prior
            LaplacesDemon::rnorminvwishart(
                mu0 = hyp$hyp$mu0[X, ],
                lambda = hyp$kappa[X],
                nu = hyp$kappa[X],
                S = hyp$hyp$Psi0[, , X] * kappa[X]
            )
        })
}

# Returns log-likelihood of single component of mixture evaluated at certain parameters
ll_normal <- function(dat, mu, sigma, prop) {
    sum(mvtnorm::dmvnorm(
        dat,
        mean = mu,
        sigma = sigma,
        log = TRUE
        ))
}

# Returns the log likelihood of the class indicators given the parameters
ll_mult <- function(dat, NIW, z, prop) {
    d <- sapply(seq_along(NIW),
                function(X)
                    mvtnorm::dmvnorm(dat,
                                     mean = NIW[[X]]$mu,
                                     sigma = NIW[[X]]$Sigma) * prop[X])
    probs <- apply(d, 1, function(X) {
        if (all(X == 0)) {
            rep(1, length(X))
        } else {
            X
        }
    })
    zeros <- rep(0, ncol(d))
    sapply(seq_along(z), function(X)
        dmultinom(replace(zeros, z[X], 1), prob = probs[,X], log = TRUE))
}

# Evaluates the density of current prop estimate at the prior
log_prop_prior <- function(prop, alpha) {
    nimble::ddirch(prop, alpha, log = TRUE)
}

# Returns log(prior * likelihood)
log_posterior <-
    function(dat,
             z,
             prop,
             NIW,
             alpha,
             mu0,
             kappa0,
             nu0,
             Psi0) {
        lprop <- log_prop_prior(prop, alpha)
        ll_normal <-
            sapply(seq_along(NIW), function(X)
                if(sum(z == X) >= ncol(dat)) {
                    ll_normal(dat[z == X,], NIW[[X]]$mu, NIW[[X]]$Sigma)    
                } else {
                    0
                })
        ll_mult <- ll_mult(dat, NIW, z, prop)
        
        lprop + sum(ll_normal) + sum(ll_mult)
    }

# Compute hastings ratio for dirichlet and NIW proposals
get_lhastings <-
    function(alpha_old,
             alpha_new) {
        alpha_old_given_new <-
            nimble::ddirch(alpha_old, alpha_new, log = TRUE)
        alpha_new_given_old <-
            nimble::ddirch(alpha_new, alpha_old, log = TRUE)
        
        alpha_old_given_new - alpha_new_given_old
    }

run_MCMC <-
    function(dat,
             init_prop,
             init_NIW,
             init_z,
             alpha,
             mu0,
             kappa0,
             Psi0,
             opt_rt = 0.3,
             iterations = 1000,
             stabilizer = 50) {
        
        # Put hyperparameters in form to match old gibbs sampler
        hyp <- list("kappa" = kappa0, "hyp" = list("mu0" = mu0, "Psi0" = Psi0))
        chain <- list()
        chain[[1]] <- list(
            "prop" = init_prop,
            "NIW" = init_NIW,
            "z" = init_z,
            "tune_var" = 2,
            "acpt" = 1
        )
        
        c0 <- 10
        c1 <- 0.8
        tune_k <- 2
        win_len <- min(iterations, 50)
        acpt <- c(1, rep(NA, win_len - 1))
        
        pb <- txtProgressBar(min = 2,
                             max = iterations,
                             style = 3)
        for (i in 2:iterations) {
            gamma1 <- c0 / (i + tune_k) ^ (c1)
            tune_var <- update_var(
                chain[[i-1]]$tune_var,
                acpt_rt = mean(acpt, na.rm = TRUE),
                opt_rt = opt_rt,
                gamma1 = gamma1
            )
            
            NIW <- draw_NIW2(dat, hyp, chain[[i - 1]]$z)
            z <- updatez(dat, NIW, chain[[i - 1]]$prop)
            prop <- propose_prop(chain[[i - 1]]$prop, tune_var = tune_var)
            
            lpost_oldprop <-
                log_posterior(
                    dat = dat,
                    z = z,
                    prop = chain[[i - 1]]$prop,
                    NIW = NIW,
                    alpha = alpha,
                    mu0 = mu0,
                    kappa0 = kappa0,
                    nu0 = nu0,
                    Psi0 = Psi0
                )
            
            lpost_newprop <-
                log_posterior(
                    dat = dat,
                    z = z,
                    prop = prop,
                    NIW = NIW,
                    alpha = alpha,
                    mu0 = mu0,
                    kappa0 = kappa0,
                    nu0 = nu0,
                    Psi0 = Psi0
                )
            
            hastings <-
                get_lhastings(
                    alpha_old = chain[[i - 1]]$prop,
                    alpha_new = prop
                )
            
            accept_ratio <-
                exp(lpost_newprop - lpost_oldprop + hastings)
            
            if (runif(1) < accept_ratio) {
                acpt[i %% win_len + 1] <- 1
                chain[[i]] <- list(
                    "prop" = prop,
                    "NIW" = NIW,
                    "z" = z,
                    "tune_var" = tune_var,
                    "acpt" = mean(acpt, na.rm = TRUE)
                )
            } else {
                acpt[i %% win_len + 1] <- 0
                chain[[i]] <- list(
                    "prop" = chain[[i - 1]]$prop,
                    "NIW" = NIW,
                    "z" = z,
                    "tune_var" = tune_var,
                    "acpt" = mean(acpt, na.rm = TRUE)
                )
            }
            setTxtProgressBar(pb, i)
        }
        close(pb)
        chain
    }