# generate data for general GMCM
rGMCM <- function(n, prop, mu, sigma){
    # n: sample size
    # prop: cluster mixing proportions
    # k: number of clusters
    # d: dimension of data (number of replicates)
    # mu: list of mean vectors of clusters
    # sigma: list of variance covariance matrices of clusters
    
    prop <- prop/sum(prop)
    k <- length(prop)
    d <- ifelse(length(prop) == 1, nrow(sigma), nrow(sigma[[1]]))
    num <- round(n*prop)
    tag <- rep(1:k, times = num)
    y <- lapply(seq(k),
                function(X) rmvnorm(num[X], mean = mu[[X]], sigma = sigma[[X]])) %>%
        abind(along = 1) %>%
        data.frame %>%
        setNames(sapply(seq(d), function(X) paste0("y.", X)))
    
    # Rescale empirical marginal distribution functions to avoid infinities
    u <- apply(y, 2, function(X) rank(X)/(n+1)) %>%
        data.frame %>%
        setNames(sapply(seq(d), function(X) paste0("u.", X)))
    
    # Calculate corresponding pseudo-data
    z <- lapply(seq(d), function(X)
        get_pseudo(u[,X],
                   map(mu,X) %>% unlist,
                   map(sigma, diag) %>% map(X) %>% unlist,
                   prop, k)) %>%
        abind(along = 2) %>%
        as.data.frame %>%
        setNames(sapply(seq(d), function(X) paste0("z.", X)))
    
    list("data" = y, "u" = u, "z" = z, "cluster" = tag)
}

# generate data for constrained version of GMCM
rconstr_GMCM <- function(n, prop, mu, sigma, rho, d){
    # n: sample size
    # prop: mixing proportion of each cluster
    # d: dimension of data (number of replicates)
    # mu: mean of "reproducible" components
    # sigma: variance of "reproducible" components
    # rho: correlation between replicates of same component
    # d: number of replicates
    # k: number of clusters
    
    # Generate all combinations of replication, given in any order
    combos <- expand.grid(rep(list(-1:1),d)) %>% as.matrix
    k <- nrow(combos)

    if(k != length(prop)){
        stop("length(prop) must be equal to total number of clusters.")
    }
    
    sigma <- c(sigma,1,sigma)
    mu <- c(-mu,0,mu)
    prop <- prop/sum(prop)
    num <- round(n*prop)
    tag <- rep(1:k, times = num)
    y <- lapply(seq(k),
                function(X) rmvnorm(num[X], mean = mu[combos[X,]+2], sigma = get_constr_sigma(diag(sigma[combos[X,]+2]), rho, combos[X,]))) %>%
        abind(along = 1) %>%
        data.frame %>%
        setNames(sapply(seq(d), function(X) paste0("y.", X)))
    
    # Rescale empirical marginal distribution functions to avoid infinities
    u <- apply(y, 2, function(X) rank(X)/(n+1)) %>%
        data.frame %>%
        setNames(sapply(seq(d), function(X) paste0("u.", X)))
    
    # Calculate corresponding pseudo-data
    z <- lapply(seq(d), function(X)
        get_pseudo(u[,X],
                   mu[combos[,X]+2],
                   sigma[combos[,X]+2],
                   prop, k)) %>%
        abind(along = 2) %>%
        as.data.frame %>%
        setNames(sapply(seq(d), function(X) paste0("z.", X)))
    
    list("data" = y, "u" = u, "z" = z, "cluster" = tag)
}
