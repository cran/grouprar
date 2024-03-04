\name{DLRule}
\alias{DLRule}
\title{
Drop the loser rule
}
\description{
Simulating drop the loser rule procedure with two-sided hypothesis testing in a clinical trial context.
}
\usage{
DLRule(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05)
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
Drop the loser rule can be describe as follows: An urn contains three types of balls \eqn{(A, B, 0)} initially. Balls of types \eqn{A} and \eqn{B} represent treatments \eqn{A} and \eqn{B}, balls of 0 type are immigration balls. If \eqn{A} (or \eqn{B}) is drawn, then treatment \eqn{A} (or \eqn{B}) is assigned to the subject and the response is observed. If the observed response is a failure, then the ball is not replaced, else replaced. If an immigration ball (type 0) is drawn, no treatment is assigned, and the ball is returned to the urn together with one \eqn{A} and one \eqn{B} ball.
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
Ivanova, A. (2003). \emph{A play-the-winner-type urn design with reduced variability.} Metrika, 58, 1-13.
}

\examples{
## a simple use
dl.res = DLRule(k = 2, p = c(0.7, 0.8), ssn = 400, Y0 = NULL, nsim = 200, alpha = 0.05)

## view the output
dl.res
\donttest{

  ## view all simulation settings
  dl.res$name
  dl.res$parameter
  dl.res$k

  ## View the simulations results
  dl.res$propotion
  dl.res$failRate
  dl.res$pwCalc
  dl.res$assignment
  }
}
\keyword{DLRlue}
\keyword{Drop the loser rule}
