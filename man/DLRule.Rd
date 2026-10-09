\name{DLRule}
\alias{DLRule}
\title{
Drop the loser rule
}
\description{
Simulating the drop-the-loser rule (number of arms \eqn{\ge 2}) with two-sided hypothesis testing in a clinical trial context.
}
\usage{
DLRule(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05, test.fun = NULL,
       typeI = FALSE, seed = NULL)
}
\arguments{
  \item{k}{A positive integer. The number of treatment groups in the trial (\eqn{k \ge 2}).}
  \item{p}{A vector of length \code{k} with values between 0 and 1. The true success rates of the treatments, used to generate data for the simulations.}
  \item{ssn}{A positive integer. The total number of participants in each simulated trial.}
  \item{Y0}{A vector of length \code{k} giving the initial number of treatment balls of each type in the urn. One immigration ball is added automatically. For instance, if \code{Y0 = c(1, 1)}, the urn starts with one ball of each treatment and one immigration ball. If \code{Y0} is \code{NULL} (default), it is set to a vector of \code{k} ones.}
  \item{nsim}{A positive integer. The number of simulated trials, with a default value of 2000.}
  \item{alpha}{A number between 0 and 1. The significance level of the two-sided test, with a default value of 0.05.}
  \item{test.fun}{An optional function \code{function(outcome, assignment)} that returns the p-value of a test for one simulated trial, where \code{outcome} holds the observed responses and \code{assignment} the treatment labels (integers 1 to \code{k}); missing responses are removed first. The null hypothesis is rejected when the p-value is at most \code{alpha}. The default \code{NULL} uses the built-in two-sided tests: the t-test for two arms and the Wald chi-squared test of equal means for more than two arms (for binary responses, if some arms have no variability, the variances are computed from \eqn{(S_k + 1)/(n_k + 2)}{(S_k + 1)/(n_k + 2)} as in Agresti and Caffo, 2000).}
  \item{typeI}{Logical. If \code{TRUE}, the simulation is repeated under the null hypothesis, with every arm's success rate set to the average of \code{p}, and the rejection rate is reported as \code{type I error}. The default is \code{FALSE}.}
  \item{seed}{An optional integer passed to \code{\link{set.seed}} before simulating, so that the results can be reproduced. The default \code{NULL} leaves the random number generator unchanged.}
}
\details{
The drop-the-loser rule can be described as follows. An urn initially contains balls of \eqn{K} treatment types and an immigration ball (type 0). A ball is drawn at random. If a treatment ball of type \eqn{k} is drawn, treatment \eqn{k} is assigned to the subject and the response is observed. If the response is a failure, the ball is not replaced, otherwise it is replaced. If the immigration ball is drawn, no treatment is assigned, and the ball is returned to the urn together with one ball of each treatment type. With two treatments \eqn{A} and \eqn{B} this is the rule of Ivanova (2003).
}
\value{
An object of class \code{"grouprar"}, a list that is printed as a short summary (see \code{\link{print.grouprar}}), with the following elements.
\item{method}{The name of the procedure.}
\item{sample size}{The total sample size.}
\item{parameter}{The true success rates used in the simulations, named \code{pA}, \code{pB}, ...}
\item{propotion}{The mean allocation proportion of each arm over the simulations, named \code{treatment A}, \code{treatment B}, ...}
\item{sd of propotion}{The standard deviation of the allocation proportion of each arm over the simulations.}
\item{failure rate}{The mean failure rate over the simulations.}
\item{sd of failure rate}{The standard deviation of the failure rate (or mean response) over the simulations.}
\item{power}{The proportion of simulated trials that reject the null hypothesis of equal success rates. Simulations in which the test cannot be computed are dropped.}
\item{data: failureRate}{The failure rate (or mean response) of each simulated trial.}
\item{data: test}{The test decision of each simulated trial (1 = reject).}
\item{data: assignment}{The treatment assignments of the last simulated trial.}
\item{data: propotion}{A data frame with the allocation proportions of each simulated trial.}
\item{data: allocation}{An \code{nsim} by \code{ssn} matrix with the treatment assignments of every simulated trial.}
\item{type I error}{Only if \code{typeI = TRUE}. The rejection rate under the null hypothesis.}
}
\references{
Ivanova, A. (2003). \emph{A play-the-winner-type urn design with reduced variability.} Metrika, 58, 1-13. \doi{10.1007/s001840200220}
}

\examples{
## a simple use
dl.res <- DLRule(k = 2, p = c(0.7, 0.8), ssn = 200, Y0 = NULL, nsim = 100, alpha = 0.05)

## view the output
dl.res
\donttest{

  ## view all simulation settings
  dl.res[["method"]]
  dl.res[["parameter"]]

  ## view the simulation results
  dl.res[["propotion"]]
  dl.res[["failure rate"]]
  dl.res[["power"]]
  dl.res[["data: assignment"]]
  }
}
\keyword{design}
\concept{drop-the-loser rule}
\concept{response-adaptive randomization}
