\name{rconstr0_GMCM}
\alias{rconstr0_GMCM}
\title{Simulate data from a Gaussian mixture copula model with relaxed constraints on behavior across replicates}
\description{
This function will simulate data from a GMCM with n data points divided up between the number of clusters determined by the number of entries in prop. Observations from a single replicate are assumed to have mean mu and variance sigma. The data have dimension d (number of replicates), and number of components h (number of possible classifications, for a given replicate). Two replicates from the same component have correlation rho. If two replicates are from different components, their correlation is 0. In \code{\link{rconstr_GMCM}}, all clusters have the same mean and correlation, up to a sign. Here, means and correlations are not necessarily equal, but must satisfy sign constraints. Variances also are permitted to vary.
}
\usage{
rconstr_GMCM(n, prop, mu, sigma, rho, d)
}
\arguments{
    \item{n}{The number of observations.}
    \item{prop}{The proportion of data that belong to each cluster. The length of prop must equal the correct number of clusters. If prop doesn't sum to 1, it is normalized to do so.}
    \item{mu}{Numeric 2-vector corresponding to the means of the negative, positive association components. The first element must be negative and the second must be positive.}
    \item{sigma}{Numeric 2-vector corresponding to the variances of the negative, positive association components.}
    \item{rho}{Numeric 3-vector corresponding to the correlation between two negative/cross/two possitive components. The first and third elements must be positive, the second must be negative.}
    \item{d}{Number of replicates.}
}
\value{
    \item{data}{An n by d numeric matrix of simulated multivariate normal data, where d is the dimension of the data.}
    \item{u}{An n by d numeric matrix corresponding to the data transformed using the empirical cdf of the observations.}
    \item{z}{An n by d numeric matrix corresponding to the latent process.}
    \item{cluster}{Integer vector of the true cluster labels.}
}
\author{\email{hbk5086@psu.edu}}
\seealso{\link{rGMCM}, \link{rconstr_GMCM}}
\examples{
set.seed(123)
pal <- sample(get_pals(4), 9, replace = FALSE)

sim <- rconstr0_GMCM(3000, rep(1/9, 9), c(-2,3), c(.9,.9), c(.8,-.5,.7), 2)

par(mfrow = c(1,2))
plot(sim$data, col = pal[sim$cluster], main = "observed")
plot(sim$u, col = pal[sim$cluster], main = "copula")
}
