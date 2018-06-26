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

# fit the general pGMCM
fpGMCM <- function(x, kmax, lambda=NULL, tol=1e-06, stepmax=50, itermax=200){
    # x: a matrix of data with rows for observations and columns for features
    # kmax: max number of clusters
    n <- nrow(x)   # sample size
    d <- ncol(x)   # dimension
    
    # initialize with pGMM
    init <- fpGMM(x, kmax, lambda = c(.1,0,1), tol=1e-04, itermax=200)
    prop0 <- init$prop
    mu0 <- init$mu
    sigma0 <- init$sigma
    k <- init$k

    # Rescale empirical marginal distribution functions to avoid infinities
    u <- apply(x, 2, function(X) rank(X)/(n+1)) %>%
        data.frame %>%
        setNames(sapply(seq(d), function(X) paste0("u.", X)))
    
    # Calculate corresponding pseudo-data
    z <- lapply(seq(d), function(X)
        get_pseudo(u[,X],
                   mu0[,X],
                   sigma0[X,X,],
                   prop0, k)) %>%
        abind(along = 2) %>%
        as.data.frame %>%
        setNames(sapply(seq(d), function(X) paste0("z.", X)))
    
    delta <- 1
    tol <- 1e-4
    ll_old <- -Inf
    
    for(step in seq(stepmax)){
        # estimation and model selection of penalized GMM for optimal lambda
        temp_fit <- fpGMM(z, kmax, lambda = NULL)
        
        # make updates
        k <- temp_fit$k
        zprop <- temp_fit$prop
        zmu <- temp_fit$mu
        zsigma <- temp_fit$sigma
        tag <- temp_fit$cluster
        zlambda <- temp_fit$lambda
        ll_new <- temp_fit$ll
        
        if(k <= 1) { break }
        
        # update pseudo-data
        z <- lapply(seq(d), function(X)
            get_pseudo(u[,X],
                       zmu[,X],
                       zsigma[X,X,],
                       zprop, k)) %>%
            abind(along = 2) %>%
            as.data.frame %>%
            setNames(sapply(seq(d), function(X) paste0("z.", X)))
        
        # measure the difference between two iteration
        delta <- abs((ll_new - ll_old) / ll_old)
        if(is.na(delta)){ delta <- 1 }

        ll_old <- ll_new
        
        if(delta < tol){ break }
        if(step > stepmax){ break }
    }

    xm <- by(data.frame(x, tag), INDICES = tag, FUN = colMeans) %>%
        abind(along = 2) %>%
        t %>%
        data.frame %>%
        mutate(tag = NULL) %>%
        unname
    
    xv <- lapply(sort(unique(tag)), function(X) var(x[tag == X,])) %>% simplify2array
    
    list("k" = k, "prop" = zprop, "mu" = xm, "sigma" = xv, "cluster" = tag, "lambda" = lambda, "ll" = ll_old)
}
