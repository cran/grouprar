\name{sqMonitor}
\alias{sqMonitor}
\alias{sqBoundary}
\title{
Group Sequential Monitoring of Response-Adaptive Randomized Trials
}
\description{
\code{sqMonitor} specifies the interim looks and the alpha spending function of a group sequential design, to be passed as the \code{monitor} argument of \code{\link{DBCD_Bin}}, \code{\link{DBCD_Cont}}, \code{\link{Group.DBCD_Bin}} and \code{\link{Group.DBCD_Cont}}. \code{sqBoundary} computes the corresponding two-sided boundaries for the standardized test statistic.
}
\usage{
sqMonitor(t, spend = "OBF")

sqBoundary(t, alpha = 0.05, spend = "OBF")
}
\arguments{
  \item{t}{A numeric vector of information times in (0, 1], the fraction of the total sample size at which the looks are taken. The final analysis at \eqn{t = 1} is added if it is missing.}
  \item{spend}{The alpha spending function. One of \code{"OBF"} (O'Brien-Fleming-type, default), \code{"Pocock"} (Pocock-type) or \code{"Linear"}.}
  \item{alpha}{A number between 0 and 1. The overall two-sided significance level, with a default value of 0.05.}
}
\details{
The spending functions of Lan and DeMets (1983) are used, with each tail spending the one-sided function at level \eqn{\alpha/2}{alpha/2} (Proschan, Lan and Wittes, 2006), as in Zhu and Hu (2010). For one tail at level \eqn{a = \alpha/2}{a = alpha/2} they are \eqn{2\{1 - \Phi(z_{a/2}/\sqrt{t})\}}{2(1 - Phi(z_(a/2)/sqrt(t)))} for \code{"OBF"}, \eqn{a\log\{1 + (e - 1)t\}}{a log(1 + (e - 1) t)} for \code{"Pocock"} and \eqn{a t} for \code{"Linear"}. The boundary \eqn{c_j} of look \eqn{j} is chosen so that, under the null hypothesis, the probability that \eqn{|Z|} first reaches the boundary at look \eqn{j} equals the alpha spent between looks \eqn{j-1} and \eqn{j}, where \eqn{Z} at the information times \eqn{t_1 < t_2 < \cdots}{t_1 < t_2 < ...} is a standard normal process with correlation \eqn{\sqrt{t_i/t_j}}{sqrt(t_i/t_j)}. Zhu and Hu (2010) showed that the sequential test statistics of a two-arm trial randomized by the doubly adaptive biased coin design have asymptotically this structure, so these boundaries keep the type I error. The boundaries are computed by recursive numerical integration of the density of \eqn{Z\sqrt{t}}{Z sqrt(t)} over the continuation region (Armitage, McPherson and Rowe, 1969), which is fast for any number of looks.

Sequential monitoring is available for two arms only, because the theory of Zhu and Hu (2010) covers the comparison of two treatments.
}
\value{
\code{sqMonitor} returns an object of class \code{"sqMonitor"}, a list with the sorted information times \code{t} and the spending function \code{spend}.

\code{sqBoundary} returns a numeric vector with the boundary for \eqn{|Z|} at each look, named \code{look 1}, \code{look 2}, ... A boundary is \code{Inf} if no alpha is spent at that look.
}
\references{
Lan, K. K. G. and DeMets, D. L. (1983). \emph{Discrete sequential boundaries for clinical trials.} Biometrika, 70(3), 659-663. \doi{10.1093/biomet/70.3.659}

Armitage, P., McPherson, C. K. and Rowe, B. C. (1969). \emph{Repeated significance tests on accumulating data.} Journal of the Royal Statistical Society, Series A, 132(2), 235-244. \doi{10.2307/2343787}

Proschan, M. A., Lan, K. K. G. and Wittes, J. T. (2006). \emph{Statistical Monitoring of Clinical Trials: A Unified Approach.} Springer.

Zhu, H. and Hu, F. (2010). \emph{Sequential monitoring of response-adaptive randomized clinical trials.} The Annals of Statistics, 38(4), 2218-2241. \doi{10.1214/10-AOS796}
}
\seealso{
\code{\link{DBCD_Bin}}, \code{\link{DBCD_Cont}}, \code{\link{Group.DBCD_Bin}}, \code{\link{Group.DBCD_Cont}}.
}
\examples{
## the boundaries of Zhu and Hu (2010): 4.877, 2.963 and 1.969
sqBoundary(c(0.2, 0.5, 1), alpha = 0.05, spend = "OBF")
sqBoundary(c(0.2, 0.5, 1), alpha = 0.05, spend = "Pocock")

## two interim looks after 1/3 and 2/3 of the patients
m <- sqMonitor(c(1/3, 2/3), spend = "OBF")
res <- DBCD_Bin(n0 = 20, p = c(0.6, 0.8), k = 2, ssn = 150, target.alloc = "RSIHR",
                nsim = 30, monitor = m, seed = 1)
res[["stopping probability"]]
}
\keyword{design}
\concept{group sequential design}
\concept{alpha spending function}
\concept{response-adaptive randomization}
