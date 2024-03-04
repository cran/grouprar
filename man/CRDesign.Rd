\name{CRDesign}
\alias{CRDesign}
\title{Complete Randomization}
\description{
Simulating complete randomization with two-sided hypothesis testing in a clinical trial context.
}
\usage{
CRDesign(k, p, ssn, nsim = 2000, alpha = 0.05)
}
\arguments{
  \item{k}{a positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})}
  \item{p}{a positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.}
  \item{ssn}{a positive integer. The value specifies the total number of participants involved in each round of the simulation.}
  \item{nsim}{a positive integer. The value specifies the total number of simulations, with a default value of 2000.}
  \item{alpha}{An integer between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.}
}
\details{
Complete randomization: Allocating participants or subjects to different treatment groups in a clinical trial in such a way that each participant has an equal and independent chance of being assigned to any of the treatment groups.
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
\examples{
## a simple use
CR.res = CRDesign(k=3, p = c(0.7, 0.8, 0.6), nsim = 500, ssn = 400)
## view the output
CR.res
\donttest{
  ## view all simulation settings
  CR.res$name
  CR.res$parameter
  CR.res$k
  ## View the simulations results
  CR.res$propotion
  CR.res$failRate
  CR.res$pwCalc
  CR.res$assignment
  }
}
\keyword{CRDesign}
\keyword{Complete randomization}
