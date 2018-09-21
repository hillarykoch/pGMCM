\name{fconstr_pGMCM}
\alias{fconstr_pGMCM}
\title{
Fit constrained penalized Gaussian copula mixture model
}
\description{
For a given maximum number of clusters in the data, find the optimal penalization parameter lambda and the optimal clustering. Lambda penalizes the mixing proportions. Model optimality is determined by best BIC. The constraints here are that, for any number of replicates, there are 3 components h = {-1,0,1}. h=0 implies no associations, while 1 and -1 imply positive and negative associations. The mean of the no association component is restricted to 0, while the positive and negative associations have mean mu and -mu. The standard deviation of the no association component is restricted to 1, while the positive and negative associations have standard deviation sigma. For two replicates where both have positive association or both have negative association, they have correlation rho. For two replicates, one with positive association and one with negative association, they have correlation -rho. Otherwise, correlation between replicates is restricted to 0.
}
\usage{
fconstr_pGMCM(x, lambda=NULL, tol=1e-06, stepmax=50, itermax=200,
                convCrit = c("GMCM", "GMM"), penaltyType = c("SCAD", "LASSO"),
                trace_params = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An n by d numeric matrix of n observations with dimension d.}
  \item{lambda}{Penalty parameter. If unspecified, a grid is automatically generated.}
  \item{tol}{Tolerance for the stopping rule for the EM algorithm. A lower tolerance will require more iterations. Defaults to 1e-06.}
  \item{stepmax}{Maximum number of updates from the copula process. Defaults to 50.}
  \item{itermax}{Maximum number of iterations of the EM algorithm to perform in the pGMM on the latent process. Defaults to 200.}
  \item{convCrit}{Character string specifying which convergence criterion to use in the pseudo-EM algorithm. If "GMCM", the algorithm stops when the difference in copula likelihoods between two consecutive steps is less than 1e-4. If "GMM", the algorithm stops when the difference between the maximum GMM likelihood betweem two consecutive steps is less than 1e-4. Note that in both cases, the maximum likelihood fit is not necessarily the final step.}
  \item{penaltyType}{Character string specifying which term to use to penalize the cluster mixing proportions. Defaults to SCAD.}
  \item{trace_params}{Logical. Should a trace of parameters, output from each best pGMM at each step, be output? Defaults to false.}
}
\details{The model is fit using a modified expectation-maximization algorithm. Clusters are initialized using a restricted GMM with 3^d clusters. Initial values of mixing proportions, means, and variance (covariance) matrices for EM are computed from these clusters.}
\value{
\item{mu}{Cluster means.}
\item{sigma}{Cluster variance (covariance) matrices}
\item{rho}{Correlation between replicates with association.}
\item{prop}{Mixing proportions in each cluster.}
\item{k}{Optimal number of clusters.}
\item{df}{Number of free parameters in final fit.}
\item{lambda}{Lambda of optimal fit.}
\item{cluster}{Cluster labels for each of the n observations.}
\item{post_prob}{Posterior probability of each observation being classified into a given group according to the final fit.}
\item{ll_gmcm}{Log copula likelihood of the data for the optimal fit.}
}
\seealso{\link{fconstr_pGMM}, \link{fpGMCM}, \link{fconstr0_pGMCM}}
\author{\email{hbk5086@psu.edu}}
\examples{
set.seed(123)
pal <- sample(get_pals(4), 9, replace = FALSE)

sim <- rconstr_GMCM(900, prop = rep(1/9,9),
                    mu = 4, sigma = 0.9, rho = 0.6, d = 2)

################################################################################
# Not run:
# # Fit with 2 different convergence criteria
# fit1 <- fconstr_pGMCM(sim$data, tol = 1e-03, stepmax = 25, itermax = 50, convCrit = "GMCM")
# fit2 <- fconstr_pGMCM(sim$data, tol = 1e-03, stepmax = 25, itermax = 50, convCrit = "GMM")
#
#par(mfrow = c(1,3))
# plot(sim$data, col = pal[sim$cluster], pch = 16, cex = 0.75,
#     xlab = expression(x[1]), ylab = expression(x[2]), main = "true clusters")
# plot(sim$data, col = pal[fit1$cluster], pch = 16, cex = 1,
#     xlab = expression(x[1]), ylab = expression(x[2]), main = "fitted clusters 1")
# plot(sim$data, col = pal[fit2$cluster], pch = 16, cex = 1,
#     xlab = expression(x[1]), ylab = expression(x[2]), main = "fitted clusters 2")
################################################################################
}
