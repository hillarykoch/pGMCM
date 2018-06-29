\name{rGMCM}
\alias{rGMCM}
\title{Simulate data from a Gaussian mixture copula model}
\usage{
rGMCM(n, prop, mu, sigma)
}
\description{
This function will simulate data from a GMCM with n data points divided up between the number of clusters determined by the number of entries in prop. The data have mean vector mu and variance-covariance matrix sigma. The dimension of the data is determined by the length of mu (or, equivalently, the dimension of sigma).
}
\arguments{
    \item{n}{The number of observations.}
    \item{prop}{The proportion of data that belong to each cluster. The length of prop specifies the number of clusters. If prop doesn't sum to 1, it is normalized to do so.}
    \item{mu}{A list of mean vectors. Each element of the list should have the same length, where the length corresponse to the dimension of the data.}
    \item{sigma}{A list of variance covariance matrices. Each element of the list should be a square matrix of the same dimension, where the dimension of the matrix should match the length of the elements in mu.}
}
\value{
    \item{data}{An n by d numeric matrix of simulated multivariate normal data, where d is the dimension of the data.}
    \item{u}{An n by d numeric matrix corresponding to the data transformed using the empirical cdf of the observations.}
    \item{z}{An n by d numeric matrix corresponding to the latent process.}
    \item{cluster}{Integer vector of the true cluster labels.}
}
\seealso{\link{rconstr_GMCM}}
\author{\email{hbk5086@psu.edu}}
\examples{
sim <- rGMCM(1500, prop = c(1/3,1/3,1/3),
            mu = list(c(-1,1), c(1,1), c(0,-sqrt(2))),
            sigma = list(rbind(c(0.65, 0.7794), c(0.7794, 1.55)),
                        rbind(c(0.65, -0.7794), c(-0.7794, 1.55)),
                        diag(c(2, 0.2))))
str(sim)
}
