# calculate pseudo variables for general GMCM
get_pseudo <- function(u, mu, sigma, prop, m) {
    # u: a vector of rescaled empirical marginal distribution function value
    # mu: a vector of means for a given replicate, 1 for each component
    # sigma: a vector of standard deviations for a given replicate, 1 for each component
    # prop: proportion of each component in mixture model

    # Compute cdf of gaussian mixture marginally for given replicate
    m <- length(prop)
    N <- 1500 # precision of grid for linear interpolation
    grid <- seq(min(-4, -mu-4*sigma), max(3, mu+4*sigma), length.out = N)
    G <- Reduce('+', lapply(seq(m), function(X) prop[X]*pnorm(grid, mu[X], sigma[X])))

    # Do linear interpolation
    interp <- approx(x = G, y = grid, xout = u)

    interp$y
}

