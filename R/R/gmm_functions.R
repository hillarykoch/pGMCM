# K-means for initialization of modified GMM
GMM_kmeans <- function(x, k, iter.max = 30) {
    # x: a matrix of data with rows for observations and columns for features
    # k: number of clusters

    n <- nrow(x)
    d <- ncol(x)
    cl <- kmeans(x, centers = k, iter.max = iter.max)
    mu <- cl$centers
    prop <- cl$size / n

    if (d == 1) {
        sigma <- lapply(seq(k), function(X)
            var(x[cl$cluster == X])) %>%
            simplify2array
    } else{
        sigma <- lapply(seq(k), function(X)
            var(x[cl$cluster == X,])) %>%
            simplify2array
    }

    list(
        "prop" = prop,
        "mu" = mu,
        "sigma" = sigma,
        "cluster" = cl$cluster
    )
}


# choose the corresponding lambda of max BIC in penalized GMM and get best estimates
fpGMM <- function(x,
                  kmax,
                  lambda = NULL,
                  tol = 1e-06,
                  itermax = 200) {
    # x: a matrix of data with rows for observations and columns for features
    # kmax: max number of clusters
    # lambda: a parameter of penalty term

    if (is.null(lambda)) {
        # Why was this hard coded like this?
        lambda <-
            sqrt(log(nrow(x))) * 10 ^ seq(-1 , 0.5, length.out = 20) # set the range of lambda
    }

    n <- nrow(x)
    d <- ncol(x)
    df <- 1 + d + d * (d + 1) / 2

  # K-means Initialization
  init <- GMM_kmeans(x, kmax)
  prop0 <- init$prop
  mu0 <- init$mu
  sigma0 <- init$sigma

    if (sum(prop0 == 0) > 0) {
        idx <- prop0 > 0
        k <- sum(idx)
        prop0 <- prop0[idx]
        mu0 <- mu0[idx,]
        sigma0 <- sigma0[, , idx]
    }

  bestBIC <- -Inf

    # Rcpp will call x a "List" if it is a data frame
    if (is.data.frame(x)) {
        x <- as.matrix(x)
    }

    for (i in seq_along(lambda)) {
        # estimate penalized GMM for a given lambda
        curGMM <- cfpGMM(
            x = x,
            prop = prop0,
            mu = mu0,
            sigma = sigma0,
            k = kmax,
            df = df,
            lambda = lambda[i],
            citermax = itermax,
            tol = tol
        )

    # parameter estimation output
    k_temp <- curGMM$k
    ll_temp <- curGMM$ll

        BIC  <- sum(ll_temp) - k_temp * df * log(n) / 2

        # update parameters
        if (bestBIC < BIC) {
            k_out <- k_temp
            prop_out <- as.vector(curGMM$prop)
            mu_out <- curGMM$mu
            sigma_out <- curGMM$sigma
            bestBIC <- BIC
            cl_out <- as.vector(curGMM$cluster)
            bestlam <- lambda[i]
            ll_out <- curGMM$ll
            post_prob <- curGMM$post_prob
        }
    }

    list(
        "k" = k_out,
        "prop" = prop_out,
        "mu" = mu_out,
        "sigma" = sigma_out,
        "cluster" = cl_out,
        "BIC" = bestBIC,
        "lambda" = bestlam,
        "ll" = ll_out,
        "post_prob" = post_prob
    )
}


# choose the corresponding lambda of max BIC in constrained penalized GMM and get best estimates
fconstr_pGMM <-
    function(x,
             lambda = NULL,
             tol = 1e-06,
             itermax = 200,
             penaltyType = c("SCAD", "LASSO"),
             bound = 0) {
        # x: a matrix of data with rows for observations and columns for features
        # kmax: max number of clusters
        # lambda: a parameter of penalty term
        # d: number of replicates
        # h: number of components (currently, make only 3)
        if (is.null(lambda)) {
            lambda <-
                sqrt(log(nrow(x))) * 10 ^ seq(-1 , 0.5, length.out = 20) # set the range of lambda
        }

        n <- nrow(x)
        d <- ncol(x)

        # Assuming there are 3 components {-1,0,1}, then max clusters should be 3^d
        # where d is the number of replicates
        kmax <- 3 ^ d
        combos <- rep(list(-1:1), d) %>% expand.grid %>% as.matrix

        # 1, 3, or 4 df depending on whether or not we need to estimate
        # rho, sigma, and mu in this constrained setting
        df <- rep(NA, kmax)
        for (i in seq(kmax)) {
            if (length(unique(abs(combos[i, ]))) < d & any(combos[i, ] != 0)) {
                df[i] <- 4
            } else if (any(combos[i, ] != 0)) {
                df[i] <- 3
            } else{
                df[i] <- 1
            }
        }

        # K-means Initialization
        # This should actually underestimate true rho, mu and overestimate true sigma
        init <- GMM_kmeans(x, kmax)
        prop0 <- init$prop
        mu0 <- c(init$mu %>% `[` (init$mu > 0.5),
                 init$mu %>% `[` (init$mu < -0.5) %>% abs) %>%
            mean %>%
            max(.5) %>%
            `*` (combos)
        sigma0 <- apply(init$sigma, 3, diag) %>%
            mean %>%
            `*` (combos) %>%
            abs
        sigma0[sigma0 == 0] <- 1
        rho0 <- c(init$sigma[1, 2, ] %>% `[` (init$sigma[1, 2, ] > 0),
                  init$sigma[1, 2, ] %>% `[` (init$sigma[1, 2, ] < 0) %>% abs) %>%
            mean %>%
            max(.1) %>%
            min(sigma0-.00001) %>%
            matrix

        if (sum(prop0 == 0) > 0) {
            idx <- prop0 > 0
            k <- sum(idx)
            prop0 <- prop0[idx]
            mu0 <- mu0[idx,]
            sigma0 <- sigma0[idx, ]
            df <- df[idx]
            combos <- combos[idx, ]
        }

        bestBIC <- -Inf

        # Rcpp will call x a "List" if it is a data frame
        if (is.data.frame(x)) {
            x <- as.matrix(x)
        }

        LASSO <- ifelse(all(penaltyType == "LASSO"), 1, 0)
        for (i in seq_along(lambda)) {
            # estimate penalized GMM for a given lambda
            curGMM <- cfconstr_pgmm(
                x = x,
                prop = prop0,
                mu = mu0,
                sigma = sigma0,
                rho = rho0,
                combos = combos,
                k = kmax,
                df = df,
                lambda = lambda[i],
                citermax = itermax,
                tol = tol,
                LASSO = LASSO,
                bound = bound
            )

            if (!any(names(curGMM) == "optim_err")) {
                ll_temp <- curGMM$ll
                df_temp <- curGMM$df
                BIC  <- sum(ll_temp) - sum(df_temp) * log(n) / 2
            } else{
                next
            }

            # update parameters
            if (bestBIC < BIC) {
                k_out <- curGMM$k
                df_out <- df_temp
                prop_out <- as.vector(curGMM$prop)
                mu_out <- curGMM$mu
                sigma_out <- curGMM$sigma
                rho_out <- curGMM$rho
                bestBIC <- BIC
                cl_out <- as.vector(curGMM$cluster)
                bestlam <- lambda[i]
                ll_out <- ll_temp
                post_prob <- curGMM$post_prob
                combos_out <- curGMM$combos
            }
        }

        # If we never have a valid fit, just return NA for now.
        tryCatch(
            expr = list(
                "k" = k_out,
                "prop" = prop_out,
                "mu" = mu_out,
                "sigma" = sigma_out,
                "rho" = rho_out,
                "df" = df_out,
                "cluster" = cl_out,
                "BIC" = bestBIC,
                "lambda" = bestlam,
                "ll" = ll_out,
                "post_prob" = post_prob,
                "combos" = combos_out
            ),
            error = function(err) NA)
  }

# choose the corresponding lambda of max BIC in constrained penalized GMM and get best estimates
# The constraints here are that less restrictive than those in fconstr_pGMM
# fconstr0_pGMM <-
#     function(x,
#              lambda = NULL,
#              tol = 1e-06,
#              itermax = 200) {
#         # x: a matrix of data with rows for observations and columns for features
#         # kmax: max number of clusters
#         # lambda: a parameter of penalty term
#         # d: number of replicates
#         # h: number of components (currently, make only 3)
#
#         if (is.null(lambda)) {
#             lambda <-
#                 sqrt(log(nrow(x))) * 10 ^ seq(-1 , 0.5, length.out = 20) # set the range of lambda
#         }
#
#         n <- nrow(x)
#         d <- ncol(x)
#
#         # Assuming there are 3 components {-1,0,1}, then max clusters should be 3^d
#         # where d is the number of replicates
#         kmax <- 3 ^ d
#         combos <- rep(list(-1:1), d) %>% expand.grid %>% as.matrix
#
#         # 1, 3, or 4 df depending on whether or not we need to estimate
#         # rho, sigma, and mu in this constrained setting
#         df <- rep(NA, kmax)
#         for (i in seq(kmax)) {
#             if (length(unique(abs(combos[i, ]))) < d & any(combos[i, ] != 0)) {
#                 df[i] <- 4
#             } else if (any(combos[i, ] != 0)) {
#                 df[i] <- 3
#             } else{
#                 df[i] <- 1
#             }
#         }
#
#         # K-means Initialization
#         init <- GMM_kmeans(x, kmax)
#         prop0 <- init$prop
#         mu0 <- c(init$mu %>% `[` (init$mu > 0),
#                  init$mu %>% `[` (init$mu < 0) %>% abs) %>%
#             mean %>%
#             `*` (combos)
#         sigma0 <- apply(init$sigma, 3, diag) %>%
#             mean %>%
#             `*` (combos) %>%
#             abs
#         sigma0[sigma0 == 0] <- 1
#         rho0 <-
#             c(
#                 init$sigma[1, 2, ] %>% `[` (init$sigma[1, 2, ] > 0) %>% mean,
#                 init$sigma[1, 2, ] %>% `[` (init$sigma[1, 2, ] < 0) %>% mean,
#                 init$sigma[1, 2, ] %>% `[` (init$sigma[1, 2, ] > 0) %>% mean
#             )
#
#         rho0_combos <- rep(0, kmax) %>%
#             replace(apply(combos, 1, sum) == -2, rho0[1]) %>%
#             replace(apply(combos, 1, sum) == 0 &
#                         !apply(combos, 1, function(X)
#                             any(X == 0)),
#                     rho0[2]) %>%
#             replace(apply(combos, 1, sum) == 2, rho0[3])
#         Sigma0 <- lapply(1:nrow(sigma0), function(X)
#             diag(sigma0[X, ]))
#         for (i in seq_along(Sigma0)) {
#             Sigma0[[i]][1, 2] <- Sigma0[[i]][2, 1] <- rho0_combos[i]
#         }
#
#         # Eliminate empty clusters
#         if (sum(prop0 == 0) > 0) {
#             idx <- prop0 > 0
#             k <- sum(idx)
#             prop0 <- prop0[idx]
#             mu0 <- mu0[idx,]
#             sigma0 <- sigma0[idx]
#             df <- df[idx]
#         }
#
#         bestBIC <- -Inf
#
#         # Rcpp will call x a "List" if it is a data frame
#         if (is.data.frame(x)) {
#             x <- as.matrix(x)
#         }
#
#         for (i in seq_along(lambda)) {
#             # estimate penalized GMM for a given lambda
#             curGMM <- cfconstr0_pGMM(
#                 x = x,
#                 prop = prop0,
#                 mu = mu0,
#                 Sigma = simplify2array(Sigma0),
#                 combos = combos,
#                 k = kmax,
#                 df = df,
#                 lambda = lambda[i],
#                 citermax = itermax,
#                 tol = tol
#             )
#
#             ll_temp <- curGMM$ll
#             df_temp <- curGMM$df
#             BIC  <- sum(ll_temp) - sum(df_temp) * log(n) / 2
#
#             # update parameters
#             if (bestBIC < BIC) {
#                 k_out <- curGMM$k
#                 df_out <- df_temp
#                 prop_out <- as.vector(curGMM$prop)
#                 mu_out <- curGMM$mu
#                 Sigma_out <- curGMM$Sigma
#                 bestBIC <- BIC
#                 cl_out <- as.vector(curGMM$cluster)
#                 bestlam <- lambda[i]
#                 ll_out <- ll_temp
#                 post_prob <- curGMM$post_prob
#             }
#         }
#
#         list(
#             "k" = k_out,
#             "prop" = prop_out,
#             "mu" = mu_out,
#             "Sigma" = Sigma_out,
#             "df" = df_out,
#             "cluster" = cl_out,
#             "BIC" = bestBIC,
#             "lambda" = bestlam,
#             "ll" = ll_out,
#             "post_prob" = post_prob
#         )
#   }
