test_mahalanobis <- function(){
    n <- 1500
    prop <- c(1/3, 1/3, 1/3)
    mu <- list(c(1,1), c(-1,-1), c(1-sqrt(2), 3+sqrt(2)))
    sigma <- list(matrix(c(1,.75,.75,1), nrow = 2, byrow = T),
                  matrix(c(2,-.8,-.8,1), nrow = 2, byrow = T),
                  diag(1:2))

    sim <- rGMCM(n, prop, mu, sigma)$data
    checkEqualsNumeric(mahalanobis(sim, mu[[1]], sigma[[1]]), Mahalanobis(as.matrix(sim), mu[[1]], sigma[[1]]))
}

test_dmvnorm <- function() {
    n <- 1500
    prop <- c(1/3, 1/3, 1/3)
    mu <- list(c(1,1), c(-1,-1), c(1-sqrt(2), 3+sqrt(2)))
    sigma <- list(matrix(c(1,.75,.75,1), nrow = 2, byrow = T),
                  matrix(c(2,-.8,-.8,1), nrow = 2, byrow = T),
                  diag(1:2))

    sim <- rGMCM(n, prop, mu, sigma)$data

    checkEqualsNumeric(dmvnorm(sim, mu[[1]], sigma[[1]]), cdmvnorm(as.matrix(sim), mu[[1]], sigma[[1]]))
}

# R code for E step
R_Estep_pGMM <- function(dat, prop, mu, sigma){
    num <- lapply(seq(mu), function(X) dmvnorm(dat, mu[[X]], sigma[[X]]) *prop[X]) %>%
        abind(along = 2)
    denom <- rowSums(num)
    h_est <- num/denom
}
