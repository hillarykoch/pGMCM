# Transform data based on quantile function of standard normal
prep_dat <- function(data, n, type = c("PIT", "scale")) {
    if(type == "scale") {
        apply(data, 2, scale)    
    } else {
        apply(data, 2, function(X) qnorm(rank(X)/(n+1)))    
    }
}