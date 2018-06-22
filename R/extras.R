# K-means for initialization of modified GMM
GMM_kmeans <- function(x, k, iter.max = 30){
    # x: a matrix of data with rows for observations and columns for features
    # k: number of clusters

    n <- nrow(x)
    d <- ncol(x)
    cl <- kmeans(x, centers = k, iter.max = iter.max)
    mu <- cl$centers
    prop <- cl$size/n

    if(d == 1){
        sigma <- lapply(seq(k), function(X) var(x[cl$cluster == X])) %>%
            simplify2array
    } else{
        sigma <- lapply(seq(k), function(X) var(x[cl$cluster == X, ])) %>%
            simplify2array
    }

    list("prop" = prop, "mu" = mu, "sigma" = sigma, "cluster" = cl$cluster)
}


# choose the corresponding lambda of max BIC in penalized GMM and get best estimates
fpGMM <- function(x, kmax, lambda=NULL, tol = 1e-06, itermax = 300){
    # x: a matrix of data with rows for observations and columns for features
    # kmax: max number of clusters
    # lambda: a parameter of penalty term

    if(is.null(lambda)){
        # Why was this hard coded like this?
        lambda <- sqrt(log(nrow(x))) * 10 ^ seq(-1 , 0.5, length.out = 20) # set the range of lambda
    }

    n <- nrow(x)
    d <- ncol(x)
    df <- 1+d+d*(d+1)/2

    # K-means Initialization
    init <- GMM_kmeans(x, kmax)
    prop0 <- init$prop
    mu0 <- init$mu
    sigma0 <- init$sigma

    if (sum(prop0 == 0) > 0) {
        idx <- prop0 > 0
        k <- sum(idx)
        prop0 <- prop0[idx]
        mu0 <- mu0[idx, ]
        sigma0 <- sigma0[,,idx]
    }

    bestBIC <- -Inf

    # Rcpp will call x a "List" if it is a data frame
    if(is.data.frame(x)){
        x <- as.matrix(x)
    }

    for(i in seq_along(lambda)){
        # estimate penalized GMM for a given lambda
        curGMM <- cfpGMM(x=x, prop=prop0, mu = mu0,
                         sigma = sigma0, k = kmax, df = df,
                         lambda = lambda[i],
                         citermax = itermax, tol = tol)

        # parameter estimation output
        k_temp <- curGMM$k
        pdf_est_temp <- curGMM$pdf_est
        prop_temp <- as.vector(curGMM$prop)

        BIC  <- sum(curGMM$ll)-k_temp*df*log(n)/2

        # update parameters
        if (bestBIC < BIC){
            k <- k_temp
            prop <- prop_temp
            mu <- curGMM$mu
            sigma <- curGMM$sigma
            bestBIC <- BIC
            cl <- as.vector(curGMM$cluster)
            bestlam <- lambda[i]
            ll <- curGMM$ll
        }
    }

    list("k" = k, "prop" = prop, "mu" = mu, "sigma" = sigma,
         "cluster" = cl, "BIC" = bestBIC, "lambda" = bestlam, "ll" = ll)
}


# choose the corresponding lambda of max BIC in constrained penalized GMM and get best estimates
# The constraints here are that
fconstr_pGMM <- function(x, kmax=NULL, lambda=NULL, tol = 1e-06, itermax = 200){
    # x: a matrix of data with rows for observations and columns for features
    # kmax: max number of clusters
    # lambda: a parameter of penalty term
    # d: number of replicates
    # h: number of components (currently, make only 3)

    if(is.null(lambda)){
        # Why was this hard coded like this?
        lambda <- sqrt(log(nrow(x))) * 10 ^ seq(-1 , 0.5, length.out = 20) # set the range of lambda
    }

    n <- nrow(x)
    d <- ncol(x)
    # Assuming there are 3 components {-1,0,1}, then max clusters should be 3^d
    # where d is the number of replicates
    kmax <- 3^d
    combos <- rep(list(-1:1), d) %>% expand.grid %>% as.matrix

    # 1, 3, or 4 df depending on whether or not we need to estimate
    # rho, sigma, and mu in this constrained setting
    df <- rep(NA, kmax)
    for(i in seq(kmax)){
        if(length(unique(abs(combos[i,]))) < d & any(combos[i,] != 0)){
            df[i] <- 4
        } else if(any(combos[i,] != 0)){
            df[i] <- 3
        } else{
            df[i] <- 1
        }
    }

    # K-means Initialization
    # This should actually underestimate true rho, mu and overestimate true sigma
    init <- GMM_kmeans(x, kmax)
    prop0 <- init$prop
    mu0 <- c(init$mu %>% `[` (init$mu > 0),
             init$mu %>% `[` (init$mu < 0) %>% abs) %>%
        mean %>%
        `*` (combos)
    sigma0 <- apply(init$sigma, 3, diag) %>%
        mean %>%
        `*` (combos) %>%
        abs
    sigma0[sigma0 == 0] <- 1
    rho0 <- c(init$sigma[1,2,] %>% `[` (init$sigma[1,2,] > 0),
              init$sigma[1,2,] %>% `[` (init$sigma[1,2,] < 0) %>% abs) %>%
        mean %>%
        matrix

    if (sum(prop0 == 0) > 0) {
        idx <- prop0 > 0
        k <- sum(idx)
        prop0 <- prop0[idx]
        mu0 <- mu0[idx, ]
        sigma0 <- sigma0[idx,]
    }

    bestBIC <- -Inf

    # Rcpp will call x a "List" if it is a data frame
    if(is.data.frame(x)){
        x <- as.matrix(x)
    }

    for(i in seq_along(lambda)){
        # estimate penalized GMM for a given lambda
        curGMM <- cfconstr_pGMM(x=x, prop=prop0, mu = mu0,
                                sigma = sigma0, rho = rho0,
                                combos = combos,
                                k = kmax, df = df,
                                lambda = lambda[i],
                                citermax = itermax, tol = tol)

        # parameter estimation output
        k_temp <- curGMM$k
        df_temp <- curGMM$df
        pdf_est_temp <- curGMM$pdf_est
        prop_temp <- as.vector(curGMM$prop)

        BIC  <- sum(curGMM$ll)-sum(df_temp)*log(n)/2

        # update parameters
        if (bestBIC < BIC){
            k <- k_temp
            prop <- prop_temp
            mu <- curGMM$mu
            sigma <- curGMM$sigma
            bestBIC <- BIC
            cl <- as.vector(curGMM$cluster)
            bestlam <- lambda[i]
            ll <- curGMM$ll
        }
    }

    list("k" = k, "prop" = prop, "mu" = mu, "sigma" = sigma, "df" = df_temp,
         "cluster" = cl, "BIC" = bestBIC, "lambda" = bestlam, "ll" = ll)
}

