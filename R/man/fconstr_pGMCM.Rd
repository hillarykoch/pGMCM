\name{fconstr_pGMCM}
\alias{fconstr_pGMCM}
\title{
Fit constrained penalized Gaussian copula mixture model
}
\description{
For a given maximum number of clusters in the data, find the optimal penalization parameter lambda and the optimal clustering. Lambda penalizes the mixing proportions. Model optimality is determined by best BIC. The constraints here are that, for any number of replicates, there are 3 components h = {-1,0,1}. h=0 implies no associations, while 1 and -1 imply positive and negative associations. The mean of the no association component is restricted to 0, while the positive and negative associations have mean mu and -mu. The standard deviation of the no association component is restricted to 1, while the positive and negative associations have standard deviation sigma. For two replicates where both have positive association or both have negative association, they have correlation rho. For two replicates, one with positive association and one with negative association, they have correlation -rho. Otherwise, correlation between replicates is restricted to 0.
}
\usage{
fconstr_pGMCM(x, kmax, lambda=NULL, tol=1e-06, stepmax=50, itermax=200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An n by d numeric matrix of n observations with dimension d.}
  \item{lambda}{Penalty parameter. If unspecified, a grid is automatically generated.}
  \item{tol}{Tolerance for the stopping rule for the EM algorithm. A lower tolerance will require more iterations. Defaults to 1e-06.}
  \item{stepmax}{Maximum number of updates from the copula process. Defaults to 50.}
  \item{itermax}{Maximum number of iterations of the EM algorithm to perform in the pGMM on the latent process. Defaults to 200.}
}
\details{The model is fit using the expectation-maximization algorithm. Clusters are initialized using kmeans with 3^d clusters. Initial values of mixing proportions, means, and variance (covariance) matrices for EM are computed from these clusters.}
\value{
\item{k}{Optimal number of clusters.}
\item{prop}{Mixing proportions in each cluster.}
\item{mu}{Cluster means.}
\item{sigma}{Cluster variance (covariance) matrices}
\item{rho}{Correlation between replicates with association.}
\item{cluster}{Cluster labels for each of the n observations.}
\item{lambda}{Lambda of optimal fit.}
\item{ll}{Log likelihood of the data given the optimal fit.}
}
\author{\email{hbk5086@psu.edu}}
\examples{
set.seed(123)
pal <- sample(get_pals(4), 9, replace = FALSE)

sim <- rconstr_GMCM(2700, prop = rep(1/9,9),
                    mu = 4, sigma = 0.7, rho = 0.6, d = 2)
fit <- fconstr_pGMCM(sim$data, kmax = 9,
                     tol = 1e-04, stepmax = 50, itermax = 150)

par(mfrow = c(1,2))
plot(sim$data, col = pal[sim$cluster], pch = 16, cex = 0.75,
    xlab = expression(x[1]), ylab = expression(x[2]), main = "true clusters")
plot(sim$data, col = pal[fit$cluster], pch = 16, cex = 0.75,
    xlab = expression(x[1]), ylab = expression(x[2]), main = "fitted clusters")
}
