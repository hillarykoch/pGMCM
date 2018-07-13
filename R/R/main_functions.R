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

    num <- ceiling(n*prop)
    rles <- rep.int(1:k, times = num) %>%
        sample(size = n, replace = F) %>%
        sort %>%
        rle
    num <- rep(0,k)
    num[rles$values] <- rles$lengths

    tag <- rep(1:k, times = num)
    y <- lapply(rles$values,
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

    if(sum(prop == 0) > 0){
        keepidx <- prop != 0
        combos <- combos[keepidx,]
        prop <- prop[keepidx]
        k <- sum(keepidx)
    }

    num <- ceiling(n*prop)
    rles <- rep.int(1:k, times = num) %>%
        sample(size = n, replace = F) %>%
        sort %>%
        rle
    num <- rep(0,k)
    num[rles$values] <- rles$lengths
    tag <- rep.int(1:k, times = num)

    y <- lapply(rles$values,
                function(X) rmvnorm(num[X],
                                    mean = mu[combos[X,]+2],
                                    sigma = get_constr_sigma(
                                        diag(sigma[combos[X,]+2]),
                                        rho,
                                        combos[X,]))) %>%
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

# generate data for second constrained version of GMCM
rconstr0_GMCM <- function(n, prop, mu, sigma, rho, d){
    # n: sample size
    # prop: mixing proportion of each cluster
    # d: dimension of data (number of replicates)
    # mu: mean of "reproducible" components
    #        a 2-vector for negative, positive association
    # sigma: variance of "reproducible" components,
    #        a 2-vector for negative, positive association
    # rho: correlation between replicates of same component
    #        a 3-vector for negative, cross, positive association
    # d: number of replicates
    # k: number of clusters

    # Generate all combinations of replication, given in any order
    combos <- expand.grid(rep(list(-1:1),d))
    k <- nrow(combos)

    if(k != length(prop)){
        stop("length(prop) must be equal to total number of clusters.")
    }
    if(length(mu) != 2){
        stop("length(mu) must equal 2.")
    }
    if(mu[1] >= 0 | mu[2] <= 0){
        stop("mu[1] must be < 0 and mu[2] must be >= 0.")
    }
    if(length(sigma) != 2){
        stop("length(sigma) must equal 2.")
    }
    if(any(sigma <= 0)){
        stop("elements of sigma must be positive.")
    }
    if(length(rho) != 3){
        stop("length(rho) must equal 3.")
    }
    if(rho[1] < 0 | rho[3] < 0 | rho[2] > 0){
        stop("rho[1] and rho[3] must be >= 0, rho[2] must be <= 0.")
    }

    sigma <- c(sigma[1],1,sigma[2])
    mu <- c(mu[1],0,mu[2])
    prop <- prop/sum(prop)

    mu_combos <- replace(combos, combos == -1, mu[1]) %>%
        replace(combos == 1, mu[3]) %>%
        as.matrix
    sig_combos <- replace(combos, combos == -1, sigma[1]) %>%
        replace(combos == 1, sigma[3]) %>%
        replace(combos == 0, 1) %>%
        as.matrix
    rho_combos <- rep(0,k) %>%
        replace(apply(combos, 1, sum) == -2, rho[1]) %>%
        replace(apply(combos, 1, sum) == 0 &
                    !apply(combos, 1, function(X) any(X == 0)), rho[2]) %>%
        replace(apply(combos, 1, sum) == 2, rho[3])

    if(sum(prop == 0) > 0){
        keepidx <- prop != 0
        combos <- combos[keepidx,]
        prop <- prop[keepidx]
        k <- sum(keepidx)
    }

    num <- ceiling(n*prop)
    rles <- rep.int(1:k, times = num) %>%
        sample(size = n, replace = F) %>%
        sort %>%
        rle
    num <- rep(0,k)
    num[rles$values] <- rles$lengths
    tag <- rep.int(1:k, times = num)

    y <- lapply(rles$values,
                function(X)
                    rmvnorm(num[X],
                            mean = mu_combos[X,],
                            sigma = matrix(c(sig_combos[X,1], rho_combos[X],
                                             rho_combos[X], sig_combos[X,2]),
                                           nrow = 2))) %>%
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
                   mu_combos[,X],
                   sig_combos[,X],
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

    for(stp in seq(stepmax)){
        # estimation and model selection of penalized GMM for optimal lambda
        temp_fit <- fpGMM(z, kmax, lambda = lambda)

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
        if(stp > stepmax){ break }
    }

    xm <- by(data.frame(x, tag), INDICES = tag, FUN = colMeans) %>%
        abind(along = 2) %>%
        t %>%
        data.frame %>%
        mutate(tag = NULL) %>%
        unname

    xv <- lapply(sort(unique(tag)), function(X) var(x[tag == X,])) %>% simplify2array

    list("k" = k, "prop" = zprop, "mu" = xm, "sigma" = xv, "cluster" = tag, "lambda" = zlambda, "ll" = ll_old)
}

# fit the constrained pGMCM
fconstr_pGMCM <- function(x, lambda=NULL, tol=1e-06, stepmax=50, itermax=200){
    # x: a matrix of data with rows for observations and columns for features
    n <- nrow(x)   # sample size
    d <- ncol(x)   # dimension

    # initialize with pGMM
    init <- fconstr_pGMM(x, lambda = c(.1,0,1), tol=1e-04, itermax=200)
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
                   sigma0[,X],
                   prop0, k)) %>%
        abind(along = 2) %>%
        as.data.frame %>%
        setNames(sapply(seq(d), function(X) paste0("z.", X)))

    delta <- 1
    tol <- 1e-4
    ll_old <- -Inf

    for(stp in seq(stepmax)){
        # estimation and model selection of penalized GMM for optimal lambda
        temp_fit <- fconstr_pGMM(z, lambda = lambda, tol = tol, itermax = itermax)

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
                       zsigma[,X],
                       zprop, k)) %>%
            abind(along = 2) %>%
            as.data.frame %>%
            setNames(sapply(seq(d), function(X) paste0("z.", X)))

        # measure the difference between two iteration
        delta <- abs((ll_new - ll_old) / ll_old)
        if(is.na(delta)){ delta <- 1 }

        ll_old <- ll_new

        if(delta < tol){ break }
        if(stp > stepmax){ break }
    }

    list("k" = k, "prop" = zprop, "mu" = zmu, "sigma" = zsigma, "rho" = temp_fit$rho, "cluster" = tag, "lambda" = zlambda, "ll" = ll_old)
}

# fit the second constrained version of pGMCM
fconstr0_pGMCM <- function(x, lambda=NULL, tol=1e-06, stepmax=50, itermax=200){
    # x: a matrix of data with rows for observations and columns for features
    n <- nrow(x)   # sample size
    d <- ncol(x)   # dimension

    # initialize with pGMM
    init <- fconstr0_pGMM(x, lambda = c(.1,0,1), tol=1e-04, itermax=200)
    prop0 <- init$prop
    mu0 <- init$mu
    Sigma0 <- init$Sigma
    k <- init$k

    # Rescale empirical marginal distribution functions to avoid infinities
    u <- apply(x, 2, function(X) rank(X)/(n+1)) %>%
        data.frame %>%
        setNames(sapply(seq(d), function(X) paste0("u.", X)))

    # Calculate corresponding pseudo-data
    z <- lapply(seq(d), function(X)
        get_pseudo(u[,X],
                   mu0[,X],
                   Sigma0[X,X,],
                   prop0, k)) %>%
        abind(along = 2) %>%
        as.data.frame %>%
        setNames(sapply(seq(d), function(X) paste0("z.", X)))

    delta <- 1
    tol <- 1e-4
    ll_old <- -Inf

    for(stp in seq(stepmax)){
        # estimation and model selection of penalized GMM for optimal lambda
        temp_fit <- fconstr0_pGMM(z, lambda = lambda, tol = tol, itermax = itermax)

        # make updates
        k <- temp_fit$k
        zprop <- temp_fit$prop
        zmu <- temp_fit$mu
        zSigma <- temp_fit$Sigma
        tag <- temp_fit$cluster
        zlambda <- temp_fit$lambda
        ll_new <- temp_fit$ll

        if(k <= 1) { break }

        # update pseudo-data
        z <- lapply(seq(d), function(X)
            get_pseudo(u[,X],
                       zmu[,X],
                       zSigma[X,X,],
                       zprop, k)) %>%
            abind(along = 2) %>%
            as.data.frame %>%
            setNames(sapply(seq(d), function(X) paste0("z.", X)))

        # measure the difference between two iteration
        delta <- abs((ll_new - ll_old) / ll_old)
        if(is.na(delta)){ delta <- 1 }

        ll_old <- ll_new

        if(delta < tol){ break }
        if(stp > stepmax){ break }
    }

    list("k" = k, "prop" = zprop, "mu" = zmu, "Sigma" = zSigma, "rho" = temp_fit$rho, "cluster" = tag, "lambda" = zlambda, "ll" = ll_old)
}
