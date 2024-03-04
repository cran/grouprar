\name{BirthDeathUrn}
\alias{BirthDeathUrn}
\title{
Birth and Death Urn
}
\description{
Simulating birth and death urn procedure (number of arms \eqn{\ge 2}) with two-sided hypothesis testing in a clinical trial context.
}
\usage{
BirthDeathUrn(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05)
}
\arguments{
  \item{k}{a positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k = 2})}
  \item{p}{a positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.}
  \item{ssn}{a positive integer. The value specifies the total number of participants involved in each round of the simulation.}
  \item{Y0}{A vector of length \code{k}, specifying the initial probability of allocating a patient to each group. For instance, if \code{Y0 = c(1, 1, 1)}, the initial probabilities are calculated as \code{Y0 / sum(Y0)}. When \code{Y0} is \code{NULL}, the initial urn will be set as If \code{Y0} is \code{NULL}, then \code{Y0} is set to a vector of length \code{k}, with all values equal to 1 by default.}
  \item{nsim}{a positive integer. The value specifies the total number of simulations, with a default value of 2000.}
  \item{alpha}{A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.}
}
\details{
The birth and death urn works as follows: Initially an urn contains balls of \eqn{K} types and \eqn{a} immigration balls. A ball is drawn randomly with replacement. If it is an immigration ball, one ball of each type is added to the urn, no patient is treated, and the next ball is drawn. The procedure is repeated until a type \eqn{i} ball \eqn{(i = 1, \cdots, K)} is drawn. Then the subject is assigned to treatment \eqn{i}. If a success, a type \eqn{i} ball is added in the urn; if a failure, a type \eqn{i} ball is removed. (Hu and Rosenberger (2006)). More detail could be found in paper \emph{A birth and death urn for randomized clinical trials} written by Ivanova etl (2000).
}
\value{
\item{name}{The name of procedure.}
\item{parameter}{The true parameters used to do the simulations.}
\item{assignment}{The randomization sequence.}
\item{propotion}{Average allocation porpotion for each of treatment groups.}
\item{failRate}{The proportion of individuals who do not achieve the expected outcome in each simulation, on average.}
\item{pwClac}{The probability of the study to detect a significant difference or effect if it truly exists.}
\item{k}{Number of arms involved in the trial.}
}
\references{
Hu, F., & Rosenberger, W. F. (2006). \emph{The theory of response-adaptive randomization in clinical trials.} John Wiley & Sons.

Ivanova, A., Rosenberger, W. F., Durham, S. D. and Flournoy, N. (2000). \emph{A birth and death urn for randomized clinical trials.} Sankhya B 62 104-118.
}

\examples{
## a simple use
bd.res = BirthDeathUrn(k = 3, p = c(0.6, 0.7, 0.6), ssn = 400, Y0 = NULL, nsim = 200, alpha = 0.05)

## view the output
bd.res
\donttest{
  ## view all simulation settings
  bd.res$name
  bd.res$parameter
  bd.res$k

  ## View the simulations results
  bd.res$propotion
  bd.res$failRate
  bd.res$pwCalc
  bd.res$assignment
  }
}

\keyword{BirthDeathUrn}
\keyword{birth and death urn}
