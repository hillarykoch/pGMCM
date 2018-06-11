## Code an EM in R to compare to

rfpGMM <- function(x, prop, mu, sigma, k, df, citermax = 300, lambda) {
    # mu a k by d matrix
    # sigma an [d,d,k] array
    n <- nrow(x)
    d <- ncol(x)
    delta <- 1
    tol = 10^(-6)

    prop_old <- prop
    mu_old <- mu
    sigma_old <- sigma

    step <- 0
    while(step < citermax){
        ## E step
        pdf_est <- prob0 <- matrix(rep(NA, n*k), nrow = n)
        for(i in 1:k){
            tmp_mu <- mu_old[i,]
            tmp_sigma <- sigma_old[,,i]
            pdf_est[,i] <- dmvnorm(x, tmp_mu, tmp_sigma)
            prob0[,i] <- pdf_est[,i]*prop_old[i]
        }
        h_est <- prob0/rowSums(prob0)

        # M step
        mu_new <- t(h_est) %*% x / colSums(h_est)
        sigma_new <- lapply(1:nrow(mu_new), function(X)
            t(x-t(outer(mu_new[X,], rep(1,n)))) %*% diag(h_est[,X]) %*%
                (x-t(outer(mu_new[X,], rep(1,n)))) / sum(h_est[,X])) %>%
            simplify2array
        prop_new <- (colSums(h_est) - lambda*df) / (n-k*lambda*df)
        if(any(prop_new < 1e-04)){
            zeroidx <- prop_new < 1e-04
            prop_new[zeroidx] <- 0
        }
        prop_new <- prop_new/sum(prop_new)

        # Difference between iterations
        delta <- sum(abs(prop_new - prop_old))
        if(delta < tol){break}

        # Eliminate small clusters
        if(any(prop_new == 0)){
            rmidx <- prop_new == 0
            k <- sum(!rmidx)

            prop_old <- prop_new[!rmidx]
            mu_old <- mu_new[!rmidx,]
            sigma_old <- sigma_new[,,!rmidx]
            pdf_est <- pdf_est[,!rmidx]
            prob0 <- prob0[,!rmidx]
            h_est <- h_est[,!rmidx]
            delta <- 1
        } else {
            prop_old <- prop_new
            mu_old <- mu_new
            sigma_old <- sigma_new
        }

        tag <- apply(h_est, 1, which.max)

        if(k <= 1){break}

    step <- step + 1
    }

    list("k" = k, "prop" = prop_old, "mu" = mu_old, "sigma" = sigma_old,
         "pdf_est" = pdf_est, "ll" = sum(prob0), "cluster" = tag)

}
