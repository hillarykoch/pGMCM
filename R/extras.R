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
# This matches fit.pGMM
fpGMM <- function(x, lambda=NULL, kmax){
    # x: a matrix of data with rows for observations and columns for features
    # kmax: max number of clusters
    # lambda: a parameter of penalty term

    if(is.null(lambda)){
        # Why was this hard coded like this?
        lambda <- sqrt(log(nrow(x))) * 10 ^ seq(- 1 , 0.5, length.out = 20) # set the range of lambda
    }

    n <- nrow(x)
    d <- ncol(x)
    df <- 1+d+d*(d+1)/2

    # K-means Initialization
    init <- GMM_kmeans(x, kmax)
    prop0 <- init$prop
    mu0 <- init$mu
    sigma0 <- init$sigma

    if (sum(prop0 ==0) > 0) {
        idx <- prop0 > 0
        k <- sum(idx)
        prop0 <- prop0[idx]
        mu0 <- mu0[idx, ]
        sigma0 <- sigma0[,,idx]
    }

    bestBIC <- -10^6

    for(i in seq_along(lambda)){
        curGMM <- cfpGMM(x=x, prop=prop0, mu = mu0, sigma = sigma0, k = kmax, df = df, citermax = 300, lambda = lambda[i]) # estimate penalized GMM for a given lambda

        # parameter estimation output
        k <- curGMM$k
        prop <- prop0 <- as.vector(curGMM$prop)
        ym <- curGMM$mu
        yv <- curGMM$sigma
        pdf_est <- curGMM$pdf_est

        tag <- as.vector(curGMM$cluster+1) # cluster indexing starts at 0 in C++

        prob <- outer(rep(1,n), prop0, "*")*pdf_est %>% colSums
        BIC  <- (k*df-1)*log(n) - sum(log(prob)) # calculate BIC score.

        # update parameters
        if (abs(bestBIC) > abs(BIC)){
            k <- curGMM$k
            prop <- curGMM$prop
            mu <- curGMM$mu
            sigma <- curGMM$sigma
            bestBIC <- BIC
            cl <- tag
            bestlam <- lambda[i]
        }
    }

    list("k" = k, "prop" = prop, "mu" = mu, "sigma" = sigma, "cluster" = cl, "BIC" = bestBIC, "lambda" = bestlam)
}
