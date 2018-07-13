# calculate pseudo variables for general GMCM
get_pseudo <- function(u, mu, sigma, prop, k) {
    # u: a vector of rescaled empirical marginal distribution function value
    # mu: a vector of means for a given replicate, 1 for each component
    # sigma: a vector of standard deviations for a given replicate, 1 for each component
    # prop: proportion of each component in mixture model

    # Compute cdf of gaussian mixture marginally for given replicate
    N <- 2000 # precision of grid for linear interpolation
    grid <- seq(min(-4.5, -mu-4.5*sigma), max(4.5, mu+4.5*sigma), length.out = N)
    G <- Reduce('+', lapply(seq(k), function(X) prop[X]*pnorm(grid, mu[X], sqrt(sigma[X]))))

    # Do linear interpolation
    interp <- approx(x = G, y = grid, xout = u)

    interp$y
}

# Calculate variance covariance matrix for constrained GMCM
get_constr_sigma <- function(Sigma, rho, idx){
    for(i in 1:(nrow(Sigma)-1)){
        for(j in (i+1):ncol(Sigma)){
            if(idx[i] == idx[j] & idx[i] != 0){
                Sigma[i,j] <- Sigma[j,i] <- rho
            } else if(idx[i] == -idx[j] & idx[i] != 0){
                Sigma[i,j] <- Sigma[j,i] <- -rho
            }
        }
    }
    Sigma
}
