\name{fpGMCM}
\alias{fpGMCM}
\title{
Fit penalized Gaussian copula mixture model
}
\description{
For a given maximum number of clusters in the data, find the optimal penalization parameter lambda and the optimal clustering. Lambda penalized the mixing proportions. Model optimality is determined by best BIC.
}
\usage{
fpGMCM(x, kmax, lambda=NULL, tol=1e-06, stepmax=50, itermax=200)
}
\arguments{
  \item{x}{An n by d numeric matrix of n observations with dimension d.}
  \item{kmax}{Maximum number of clusters to consider, and the number of clusters the algorithm is initialized with.}
  \item{lambda}{Penalty parameter. If unspecified, a grid is automatically generated.}
  \item{tol}{Tolerance for the stopping rule for the EM algorithm. A lower tolerance will require more iterations. Defaults to 1e-06.}
  \item{stepmax}{Maximum number of updates from the copula process. Defaults to 50.}
  \item{itermax}{Maximum number of iterations of the EM algorithm to perform in the pGMM on the latent process. Defaults to 200.}
}
\details{The model is fit using the expectation-maximization algorithm. Clusters are initialized using kmeans with kmax clusters. Initial values of mixing proportions, means, and variance (covariance) matrices for EM are computed from these clusters.}
\value{
\item{k}{Optimal number of clusters.}
\item{prop}{Mixing proportions in each cluster.}
\item{mu}{Cluster means.}
\item{sigma}{Cluster variance (covariance) matrices}
\item{cluster}{Cluster labels for each of the n observations.}
\item{lambda}{Lambda of optimal fit.}
\item{ll}{Log likelihood of the data given the optimal fit.}
}
\author{\email{hbk5086@psu.edu}}
\examples{
library(mvtnorm)
set.seed(123)
pal <- get_pals(7)

n <- 1500
k <- 3
d <- 2
kmax <- 6

prop <- c(1/3, 1/3, 1/3)
mu <- list(c(-1, 1), c(1, 1), c(0, -sqrt(2)))
sigma <- list(rbind(c(0.65, 0.7794), c(0.7794, 1.55)),
                 rbind(c(0.65, -0.7794), c(-0.7794, 1.55)),
                 diag(c(2, 0.2)))

sim <- rGMCM(n, prop, mu, sigma)

################################################################################
# Not run:
# fit <- fpGMCM(sim$data, kmax = 10)
#
# par(mfrow = c(1,2))
# plot(sim$data, col = rev(pal)[sim$cluster], pch = 16, cex = 0.75,
#    xlab = expression(x[1]), ylab = expression(x[2]), main = "true clusters")
# plot(sim$data, col = rev(pal)[fit$cluster], pch = 16, cex = 0.75,
#    xlab = expression(x[1]), ylab = expression(x[2]), main = "fitted clusters")
################################################################################
}
