\name{WeiUrn}
\alias{WeiUrn}
\title{
Randomized Play-the-winner rule with multiple arms (\eqn{k > 2})
}
\description{
Simulating randomized play-the-winner rule (multiple arms) with two-sided hypothesis testing in a clinical trial context.
}
\usage{
  WeiUrn(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05)
}
\arguments{
  \item{k}{a positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k > 2})}
  \item{p}{a positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.}
  \item{ssn}{a positive integer. The value specifies the total number of participants involved in each round of the simulation.}
  \item{Y0}{A vector of length \code{k}, specifying the initial probability of allocating a patient to each group. For instance, if \code{Y0 = c(1, 1, 1)}, the initial probabilities are calculated as \code{Y0 / sum(Y0)}. When \code{Y0} is \code{NULL}, the initial urn will be set as If \code{Y0} is \code{NULL}, then \code{Y0} is set to a vector of length \code{k}, with all values equal to 1 by default.}
  \item{nsim}{a positive integer. The value specifies the total number of simulations, with a default value of 2000.}
  \item{alpha}{A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.}
}
\details{
Wei's urn procedure is obtained by extending the randomized play the winner rule (Wei1978) from the case \eqn{k = 2} to \eqn{k > 2}. Hence, It enables to conduct multi-arm clinical trials, and offers a greater range of applications.
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
LJ Wei (1979). \emph{The generalized polya’s urn design for sequential medical trials.} The Annals of Statistics, 7(2):291–296, 19
}

\examples{
## a simple use
wei.res = WeiUrn(k = 3, p = c(0.7, 0.8, 0.7), ssn = 400, Y0 = NULL, nsim = 200, alpha = 0.05)

## view the output
wei.res
\donttest{

  ## view all simulation settings
  wei.res$name
  wei.res$parameter
  wei.res$k

  ## View the simulations results
  wei.res$propotion
  wei.res$failRate
  wei.res$pwCalc
  wei.res$assignment
  }
}

\keyword{WeiUrn}
\keyword{Randomized Play-the-winner rule with multiple arms}
