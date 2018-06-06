# generate data for general GMCM
rGMCM <- function(n, prop, mu, sigma){
    # n: sample size
    # prop: proportion of each component in mixture model
    # m: number of components (e.g. h \in {-1, 0, 1} is 3 components)
    # d: dimension of data (number of replicates)
    # mu: list of mean vectors of components
    # sigma: list of variance covariance matrices of components

    prop <- prop/sum(prop)
    m <- length(prop)
    d <- ifelse(length(prop) == 1, nrow(sigma), nrow(sigma[[1]]))
    num <- round(n*prop)
    comp <- rep(1:m, times = num)
    y <- lapply(seq(m),
                function(X) rmvnorm(num[X], mean = mu[[X]], sigma = sigma[[X]])) %>%
        abind(along = 1) %>%
        data.frame %>%
        setNames(c("y.1", "y.2"))

    # Rescale empirical marginal distribution functions to avoid infinities
    u <- apply(y, 2, function(X) rank(X)/(n+1)) %>%
        data.frame %>%
        setNames(c("u.1", "u.2"))

    # Calculate corresponding pseudo-data
    z <- lapply(seq(d), function(X)
        get_pseudo(u[,X],
                   map(mu,X) %>% unlist,
                   map(sigma, diag) %>% map(X) %>% unlist,
                   prop, m)) %>%
        abind(along = 2) %>%
        as.data.frame %>%
        setNames(c("z.1", "z.2"))

    list("data" = y, "u" = u, "z" = z, "component" = comp)
}
