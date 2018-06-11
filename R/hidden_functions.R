# calculate pseudo variables for general GMCM
get_pseudo <- function(u, mu, sigma, prop, k) {
    # u: a vector of rescaled empirical marginal distribution function value
    # mu: a vector of means for a given replicate, 1 for each component
    # sigma: a vector of standard deviations for a given replicate, 1 for each component
    # prop: proportion of each component in mixture model

    # Compute cdf of gaussian mixture marginally for given replicate
    N <- 1500 # precision of grid for linear interpolation
    grid <- seq(min(-4, -mu-4*sigma), max(3, mu+4*sigma), length.out = N)
    G <- Reduce('+', lapply(seq(k), function(X) prop[X]*pnorm(grid, mu[X], sigma[X])))

    # Do linear interpolation
    interp <- approx(x = G, y = grid, xout = u)

    interp$y
}

# Calculate variance covariance matrix for constrained GMCM
get_constr_sigma <- function(sigma, idx){
    Sigma <- diag(sigma[idx])
    for(i in 1:(nrow(Sigma)-1)){
        for(j in (i+1):ncol(Sigma)){
            if(Sigma[i,i] == Sigma[j,j]){
                Sigma[i,j] <- Sigma[j,i] <- rho[idx[i]]
            }
        }
    }
    Sigma
}
