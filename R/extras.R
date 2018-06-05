# K-means for initialization of modified GMM
GMM_kmeans <- function(y, k, iter.max = 30){
    # y: a matrix of data with rows for observations and columns for features
    # k: number of clusters
    n <- nrow(y)
    d <- ncol(y)
    cl <- kmeans(y, centers = k, iter.max = iter.max)
    mu <- cl$centers
    prop <- cl$size/n
    
    if(d == 1){
        sigma <- lapply(seq(k), function(X) var(y[cl$cluster == X]))
    } else{
        sigma <- lapply(seq(k), function(X) var(y[cl$cluster == X, ]))
    }
    
    list("prop" = prop, "mu" = mu, "sigma" = sigma, "cluster" = cl$cluster)
}
