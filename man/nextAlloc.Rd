\name{nextAlloc}
\alias{nextAlloc}
\title{
Allocation Probabilities for an Ongoing Trial
}
\description{
Computes the allocation probabilities of the next patient, or of the next group of patients, in an ongoing response-adaptive randomized trial from the assignments and responses observed so far, and optionally draws the assignments. The same estimates, targets and allocation functions as in the simulation functions (for example \code{\link{DBCD_Bin}} and \code{\link{Group.DBCD_Bin}}) are used.
}
\usage{
nextAlloc(alloc, outcome, k, response = "binary", target.alloc = NULL,
          allocation = "DBCD", r = 2, erade.alpha = 0.5, lower.bound = 0,
          theta0 = NULL, size = 1)
}
\arguments{
  \item{alloc}{An integer vector with the treatment (1 to \code{k}) of every patient enrolled so far.}
  \item{outcome}{A numeric vector of the same length as \code{alloc} with the response of each patient (0 or 1 for binary responses). Use \code{NA} for responses that are not observed yet or are missing.}
  \item{k}{A positive integer. The number of treatment groups in the trial (\eqn{k \ge 2}).}
  \item{response}{The response type, \code{"binary"} (default) or \code{"continuous"}.}
  \item{target.alloc}{Desired allocation proportion. For binary responses one of \code{"Neyman"}, \code{"RSIHR"}, \code{"RPW"}, \code{"WeisUrn"}, \code{"OptimalNeyman"} or \code{"OptimalRSIHR"}, with default \code{"RPW"}. For continuous responses one of \code{"Neyman"}, \code{"ZR"}, \code{"OptimalNeyman"} or \code{"DaOptimal"}, with default \code{"Neyman"}. See \code{\link{DBCD_Bin}} and \code{\link{DBCD_Cont}} for their definitions.}
  \item{allocation}{The allocation function, \code{"DBCD"} (default, Hu and Zhang, 2004) or \code{"ERADE"} (Hu, Zhang and He, 2009; Alkhnefr, Hu and Zhai, 2025). See \code{\link{DBCD_Bin}}.}
  \item{r}{A non-negative number. The tuning parameter of Hu and Zhang's allocation function, used when \code{allocation = "DBCD"}. The default value is 2.}
  \item{erade.alpha}{A number between 0 and 1. The degree of randomization of ERADE, used when \code{allocation = "ERADE"}. The default is 0.5.}
  \item{lower.bound}{A number between 0 and \eqn{1/k}. The smallest allowed target proportion of each arm for the optimal targets and \code{"ZR"}. The default is 0.}
  \item{theta0}{A vector of length \code{k} used to smooth the success-rate estimates of binary responses, \eqn{\hat p_k = (S_k + \theta_{0k})/(N_k + 1)}{p_k = (S_k + theta0_k)/(N_k + 1)}, where \eqn{S_k} and \eqn{N_k} are the number of observed successes and observed responses on treatment \eqn{k}. If \code{NULL} (default), all values are 0.5. Not used for continuous responses.}
  \item{size}{A non-negative integer. The number of assignments to draw with the returned probabilities, for example the size of the next group. The default is 1, and 0 draws none.}
}
\details{
The parameters of each arm are estimated from the observed responses only (\code{NA} values are ignored): smoothed success rates for binary responses, and the mean and variance (which needs at least two responses) for continuous responses. The current allocation proportions use all enrolled patients, including those whose responses are not available yet. If the target cannot be estimated, equal allocation is used as the target, and with no enrolled patients every arm has probability \eqn{1/k}. All patients of a group are assigned independently with the same probabilities, as in the group designs of Zhai, Li, Zhang and Hu (2024).
}
\value{
A list with the following elements.
\item{prob}{The allocation probabilities of the next patient for each arm, named \code{treatment A}, \code{treatment B}, ...}
\item{target}{The estimated target allocation proportions.}
\item{estimate}{The parameter estimates, named \code{pA}, \code{pB}, ... or \code{muA}, \code{sigma2A}, ...}
\item{assignment}{An integer vector of \code{size} assignments drawn with probabilities \code{prob}.}
}
\references{
Alkhnefr, N., Hu, F. and Zhai, G. (2025). \emph{Efficient randomized adaptive designs for multi-arm clinical trials.} Statistical Methods in Medical Research, 34(9), 1886-1898. \doi{10.1177/09622802251362644}

Hu, F. and Zhang, L. X. (2004). \emph{Asymptotic properties of doubly adaptive biased coin designs for multitreatment clinical trials.} The Annals of Statistics, 32(1), 268-301. \doi{10.1214/aos/1079120137}

Hu, F., Zhang, L. X. and He, X. (2009). \emph{Efficient randomized-adaptive designs.} The Annals of Statistics, 37(5A), 2543-2560. \doi{10.1214/08-AOS655}

Zhai, G., Li, Y., Zhang, L. and Hu, F. (2024). \emph{Group response-adaptive randomization with delayed and missing responses.} Statistics in Medicine, 43(27), 5047-5059. \doi{10.1002/sim.10220}
}
\seealso{
\code{\link{DBCD_Bin}}, \code{\link{DBCD_Cont}}, \code{\link{Group.DBCD_Bin}} for simulating the designs.
}
\examples{
## 30 patients enrolled in a three-arm trial, the last three responses are pending
set.seed(1)
alloc <- sample(1:3, 30, replace = TRUE)
outcome <- rbinom(30, 1, c(0.6, 0.7, 0.8)[alloc])
outcome[28:30] <- NA

## probabilities for the next patient with the DBCD and the RSIHR target
nextAlloc(alloc, outcome, k = 3, target.alloc = "RSIHR")

## assignments for the next group of 5 patients with ERADE
nextAlloc(alloc, outcome, k = 3, target.alloc = "RSIHR", allocation = "ERADE",
          size = 5)$assignment

## continuous responses
y <- rnorm(30, mean = c(13, 15, 14)[alloc], sd = 3)
nextAlloc(alloc, y, k = 3, response = "continuous", target.alloc = "Neyman")$prob
}
\keyword{design}
\concept{response-adaptive randomization}
\concept{doubly adaptive biased coin design}
\concept{efficient randomized adaptive design}
