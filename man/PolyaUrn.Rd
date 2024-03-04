\name{PolyaUrn}
\alias{PolyaUrn}
\title{
Randomized Pólya urn procedure
}
\description{
Simulating randomized Pólya urn procedure with two-sided hypothesis testing in a clinical trial context.
}
\usage{
PolyaUrn(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05)
}
\arguments{
  \item{k}{a positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})}
  \item{p}{a positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.}
  \item{ssn}{a positive integer. The value specifies the total number of participants involved in each round of the simulation.}
  \item{Y0}{A vector of length \code{k}, specifying the initial probability of allocating a patient to each group. For instance, if \code{Y0 = c(1, 1, 1)}, the initial probabilities are calculated as \code{Y0 / sum(Y0)}. When \code{Y0} is \code{NULL}, the initial urn will be set as If \code{Y0} is \code{NULL}, then \code{Y0} is set to a vector of length \code{k}, with all values equal to 1 by default.}
  \item{nsim}{a positive integer. The value specifies the total number of simulations, with a default value of 2000.}
  \item{alpha}{A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.}
}
\details{
The randomized Pólya urn (RPU) procedure can be describe as follows: An urn contains at least one ball of each treatment type (totally K treatments) initially. A ball is drawn from the urn with replacement. If a type \eqn{i} ball is drawn, \eqn{i=1, \ldots, K}, then treatment \eqn{i} is assigned to the next patient. If the response is a success, a ball of type \eqn{i} is added to the urn. Otherwise the urn remains unchanged.}

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
Durham, S. D., FlournoY, N. AND LI, W. (1998). \emph{Sequential designs for maximizing the probability of a favorable response.} Canadian Journal of Statistics, 3, 479-495.
}

\examples{
## a simple use
Polya.res = PolyaUrn(k = 3, p = c(0.6, 0.7, 0.6), ssn = 400, Y0 = NULL, nsim = 200, alpha = 0.05)

## view the output
Polya.res
\donttest{
  ## view all simulation settings
  Polya.res$name
  Polya.res$parameter
  Polya.res$k

  ## View the simulations results
  Polya.res$propotion
  Polya.res$failRate
  Polya.res$pwCalc
  Polya.res$assignment
  }
}

\keyword{PolyaUrn}
\keyword{Randomized Pólya Urn}
