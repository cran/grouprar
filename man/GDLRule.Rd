\name{GDLRule}
\alias{GDLRule}
\title{
Generalized drop-the-loser rule
}
\description{
Simulating generalized drop-the-loser rule procedure (number of arms \eqn{> 2}) with two-sided hypothesis testing in a clinical trial context.
}
\usage{
GDLRule(k, p, ssn, aK, Y0 = NULL, nsim = 2000, alpha = 0.05)
}
\arguments{
  \item{k}{a positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})}
  \item{p}{a positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.}
  \item{ssn}{a positive integer. The value specifies the total number of participants involved in each round of the simulation.}
  \item{aK}{a positive vector of length equals to \code{k}. The values specifies when the immigration ball is drawn, the number of each treatment ball added to the urn.}
  \item{Y0}{A vector of length \code{k}, specifying the initial probability of allocating a patient to each group. For instance, if \code{Y0 = c(1, 1, 1)}, the initial probabilities are calculated as \code{Y0 / sum(Y0)}. When \code{Y0} is \code{NULL}, the initial urn will be set as If \code{Y0} is \code{NULL}, then \code{Y0} is set to a vector of length \code{k}, with all values equal to 1 by default.}
  \item{nsim}{a positive integer. The value specifies the total number of simulations, with a default value of 2000.}
  \item{alpha}{A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.}
}
\details{
Consider an urn containing balls of \eqn{K+1} types. Balls of types $1, \ldots, K$ represent treatments, balls of type 0 will be called immigration balls. When the subject arrived for randomizition, a ball is drawn at random. If the ball is of type 0 (i.e, an immigration ball), no subject is treated, and the ball is returned to the urn together with \eqn{A=a_1+\cdots+a_K} additional balls, \eqn{a_k} of treatment type \eqn{k, k=1, \ldots, K}. If a treatment ball is drawn (say, of type \eqn{k}, for some \eqn{k=1, \ldots, K}) the next subject is given treatment \eqn{k} and the ball is not replaced. If the observed response of this subject is a success, then the ball is replaced, otherwise not replaced.
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
Zhang, L. X., Chan, W. S., Cheung, S. H., & Hu, F. (2007). \emph{A generalized drop-the-loser urn for clinical trials with delayed responses.} Statistica Sinica, 17(1), 387-409.
}

\examples{
## a simple use
gdl.res = GDLRule(k = 3, p = c(0.6, 0.7, 0.6),
                  ssn = 400, aK = c(1, 1, 1), Y0 = NULL, nsim = 200, alpha = 0.05)

## view the output
gdl.res
\donttest{

  ## view all simulation settings
  gdl.res$name
  gdl.res$parameter
  gdl.res$k

  ## View the simulations results
  gdl.res$propotion
  gdl.res$failRate
  gdl.res$pwCalc
  gdl.res$assignment
  }
}

\keyword{GDLRlue}
\keyword{Generalized Drop the loser rule}
