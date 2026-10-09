\name{Bai.Hu.Shen.Urn}
\alias{Bai.Hu.Shen.Urn}
\title{
Bai Hu Shen's Urn
}
\description{
Bai, Hu, and Shen (2002) proposed an adaptive design for multi-arm clinical trials. The allocation probabilities adapt to the performance of the patients already treated: a success on a treatment increases the chance that the next patient is assigned to it, and a failure moves probability to the other treatments in proportion to their estimated success rates. This function simulates the Bai, Hu, and Shen urn with two-sided hypothesis testing in a clinical trial context.
}

\usage{
Bai.Hu.Shen.Urn(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05,
                test.fun = NULL, typeI = FALSE, seed = NULL)
}

\arguments{
  \item{k}{A positive integer. The number of treatment groups in the trial (\eqn{k \ge 2}).}
  \item{p}{A vector of length \code{k} with values between 0 and 1. The true success rates of the treatments, used to generate data for the simulations.}
  \item{ssn}{A positive integer. The total number of participants in each simulated trial.}
  \item{Y0}{A vector of length \code{k} giving the initial urn composition (number of balls of each treatment type). For instance, if \code{Y0 = c(1, 1, 1)}, the first patient is assigned to each treatment with probability \code{Y0 / sum(Y0)}. If \code{Y0} is \code{NULL} (default), it is set to a vector of \code{k} ones.}
  \item{nsim}{A positive integer. The number of simulated trials, with a default value of 2000.}
  \item{alpha}{A number between 0 and 1. The significance level of the two-sided test, with a default value of 0.05.}
  \item{test.fun}{An optional function \code{function(outcome, assignment)} that returns the p-value of a test for one simulated trial, where \code{outcome} holds the observed responses and \code{assignment} the treatment labels (integers 1 to \code{k}); missing responses are removed first. The null hypothesis is rejected when the p-value is at most \code{alpha}. The default \code{NULL} uses the built-in two-sided tests: the t-test for two arms and the Wald chi-squared test of equal means for more than two arms (for binary responses, if some arms have no variability, the variances are computed from \eqn{(S_k + 1)/(n_k + 2)}{(S_k + 1)/(n_k + 2)} as in Agresti and Caffo, 2000).}
  \item{typeI}{Logical. If \code{TRUE}, the simulation is repeated under the null hypothesis, with every arm's success rate set to the average of \code{p}, and the rejection rate is reported as \code{type I error}. The default is \code{FALSE}.}
  \item{seed}{An optional integer passed to \code{\link{set.seed}} before simulating, so that the results can be reproduced. The default \code{NULL} leaves the random number generator unchanged.}
}
\details{
Bai, Hu and Shen's urn can be described as follows. An urn initially contains balls of \eqn{K} types, where balls of types \eqn{1, 2, \ldots, K} represent treatments \eqn{1, 2, \ldots, K}. A ball is drawn at random from the urn and, if it is of type \eqn{k}, the next patient is assigned to treatment \eqn{k}. After the response is observed, the urn composition is updated. A success on treatment \eqn{k} adds one ball of type \eqn{k} to the urn. A failure on treatment \eqn{k} adds \eqn{\hat p_j/(\hat M - \hat p_k)}{p_j/(M - p_k)} balls of each other type \eqn{j \ne k}{j != k}, where \eqn{\hat p_j = (S_j + 1)/(N_j + 1)}{p_j = (S_j + 1)/(N_j + 1)} is the current estimate of the success rate of treatment \eqn{j} (\eqn{S_j} successes among \eqn{N_j} patients) and \eqn{\hat M = \hat p_1 + \cdots + \hat p_K}{M = p_1 + ... + p_K}. The estimates use the responses of the previous patients only. This is adaptive design 3 of Bai, Hu and Shen (2002), the design proposed in the paper. (Versions of grouprar before 0.2.0 used the true success rates instead, which is their design 2.)
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
Bai, Z. D., Hu, F. and Shen, L. (2002). \emph{An adaptive design for multi-arm clinical trials.} Journal of Multivariate Analysis, 81(1), 1-18. \doi{10.1006/jmva.2001.1987}
}

\examples{
## a simple use
bhs.res <- Bai.Hu.Shen.Urn(k = 3,
                           p = c(0.7, 0.8, 0.6),
                           ssn = 200,
                           Y0 = NULL,
                           nsim = 100,
                           alpha = 0.05)

## view the output
bhs.res
\donttest{
  ## view all simulation settings
  bhs.res[["method"]]
  bhs.res[["parameter"]]

  ## view the simulation results
  bhs.res[["propotion"]]
  bhs.res[["failure rate"]]
  bhs.res[["power"]]
  bhs.res[["data: assignment"]]
  }
}
\keyword{design}
\concept{Bai Hu Shen urn}
\concept{response-adaptive randomization}
