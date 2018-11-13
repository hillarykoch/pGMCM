# Transform data based on quantile function of standard normal
prep_dat <- function(data, n) {
    apply(data, 2, function(X) qnorm(rank(X)/(n+1)))
}