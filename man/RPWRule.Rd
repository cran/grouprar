\name{RPWRule}
\alias{RPWRule}
\title{
Randomized Play-the-winner Rule
}
\description{
Simulating randomized play-the-winner rule with two-sided hypothesis testing in a clinical trial context.
}
\usage{
RPWRule(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05)
}
\arguments{
   \item{k}{a positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k = 2})}
  \item{p}{a positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.}
  \item{ssn}{a positive integer. The value specifies the total number of participants involved in each round of the simulation.}
  \item{Y0}{A vector of length \code{k}, specifying the initial probability of allocating a patient to each group. For instance, if \code{Y0 = c(1, 1,)}, the initial probabilities are calculated as \code{Y0 / sum(Y0)}. When \code{Y0} is \code{NULL}, the initial urn will be set as If \code{Y0} is \code{NULL}, then \code{Y0} is set to a vector of length \code{k}, with all values equal to 1 by default.}
  \item{nsim}{a positive integer. The value specifies the total number of simulations, with a default value of 2000.}
  \item{alpha}{An integer between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.}
}
\details{
The Randomized Play-the-Winner Rule allocates future subjects in a clinical trial to treatment groups based on the performance of previously treated subjects. This rule increases the likelihood of future patients being assigned to the better-performing treatment, as determined by the outcomes of previously treated subjects.
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
L. J. WEI and S. DURHAM (1978) \emph{The Randomized Play-the-Winner Rule in Medical Trials.} Journal of the American Statistical Association, 73, 364, 840--843.
}

\examples{
## a simple use
RPW.res = RPWRule(k = 2, p = c(0.7, 0.8), ssn = 400, Y0 = NULL, nsim = 200, alpha = 0.05)
## view the output
RPW.res
\donttest{

  ## view all simulation settings
  RPW.res$name
  RPW.res$parameter
  RPW.res$k

  ## View the simulations results
  RPW.res$propotion
  RPW.res$failRate
  RPW.res$pwCalc
  RPW.res$assignment
  }
}

\keyword{RPWRule}
\keyword{Randomized Play-the-winner Rule}
