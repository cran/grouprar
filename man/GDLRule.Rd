\name{GDLRule}
\alias{GDLRule}
\title{
Generalized drop-the-loser rule
}
\description{
Simulating the generalized drop-the-loser rule (number of arms \eqn{\ge 2}) with two-sided hypothesis testing in a clinical trial context.
}
\usage{
GDLRule(k, p, ssn, aK, Y0 = NULL, nsim = 2000, alpha = 0.05, test.fun = NULL,
        typeI = FALSE, seed = NULL)
}
\arguments{
  \item{k}{A positive integer. The number of treatment groups in the trial (\eqn{k \ge 2}).}
  \item{p}{A vector of length \code{k} with values between 0 and 1. The true success rates of the treatments, used to generate data for the simulations.}
  \item{ssn}{A positive integer. The total number of participants in each simulated trial.}
  \item{aK}{A positive vector of length \code{k}. The number of balls of each treatment type added to the urn when the immigration ball is drawn.}
  \item{Y0}{A vector of length \code{k} giving the initial number of treatment balls of each type in the urn. One immigration ball is added automatically. For instance, if \code{Y0 = c(1, 1, 1)}, the urn starts with one ball of each treatment and one immigration ball. If \code{Y0} is \code{NULL} (default), it is set to a vector of \code{k} ones.}
  \item{nsim}{A positive integer. The number of simulated trials, with a default value of 2000.}
  \item{alpha}{A number between 0 and 1. The significance level of the two-sided test, with a default value of 0.05.}
  \item{test.fun}{An optional function \code{function(outcome, assignment)} that returns the p-value of a test for one simulated trial, where \code{outcome} holds the observed responses and \code{assignment} the treatment labels (integers 1 to \code{k}); missing responses are removed first. The null hypothesis is rejected when the p-value is at most \code{alpha}. The default \code{NULL} uses the built-in two-sided tests: the t-test for two arms and the Wald chi-squared test of equal means for more than two arms (for binary responses, if some arms have no variability, the variances are computed from \eqn{(S_k + 1)/(n_k + 2)}{(S_k + 1)/(n_k + 2)} as in Agresti and Caffo, 2000).}
  \item{typeI}{Logical. If \code{TRUE}, the simulation is repeated under the null hypothesis, with every arm's success rate set to the average of \code{p}, and the rejection rate is reported as \code{type I error}. The default is \code{FALSE}.}
  \item{seed}{An optional integer passed to \code{\link{set.seed}} before simulating, so that the results can be reproduced. The default \code{NULL} leaves the random number generator unchanged.}
}
\details{
Consider an urn containing balls of \eqn{K+1} types. Balls of types \eqn{1, \ldots, K} represent treatments, and balls of type 0 are called immigration balls. When a subject arrives for randomization, a ball is drawn at random. If the ball is of type 0 (an immigration ball), no subject is treated, and the ball is returned to the urn together with \eqn{A=a_1+\cdots+a_K} additional balls, \eqn{a_k} of treatment type \eqn{k, k=1, \ldots, K}. If a treatment ball is drawn (say, of type \eqn{k}, for some \eqn{k=1, \ldots, K}), the subject is given treatment \eqn{k}. If the observed response of this subject is a success, the ball is replaced, otherwise it is not replaced.
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
Zhang, L. X., Chan, W. S., Cheung, S. H. and Hu, F. (2007). \emph{A generalized drop-the-loser urn for clinical trials with delayed responses.} Statistica Sinica, 17(1), 387-409.
}

\examples{
## a simple use
gdl.res <- GDLRule(k = 3, p = c(0.6, 0.7, 0.6),
                   ssn = 200, aK = c(1, 1, 1), Y0 = NULL, nsim = 100, alpha = 0.05)

## view the output
gdl.res
\donttest{

  ## view all simulation settings
  gdl.res[["method"]]
  gdl.res[["parameter"]]

  ## view the simulation results
  gdl.res[["propotion"]]
  gdl.res[["failure rate"]]
  gdl.res[["power"]]
  gdl.res[["data: assignment"]]
  }
}

\keyword{design}
\concept{generalized drop-the-loser rule}
\concept{response-adaptive randomization}
