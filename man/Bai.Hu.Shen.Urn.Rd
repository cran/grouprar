\name{Bai Hu Shen's Urn}
\alias{Bai.Hu.Shen.Urn}
\title{
Bai Hu Shen's Urn
}
\description{
Bai, Hu, and Shen \cite{bai2002adaptive} proposed a new adaptive design for multi-arm clinical trials. The main idea behind this procedure is that the allocation probability adapts based on the performance of the most recent patients under their assigned treatment. Positive performance in that treatment increase the likelihood of the next patient being assigned to this group, whereas negative outcomes decrease it. This function is for simulating the Bai, Hu, and Shen's urn model under two-sided hypothesis testing in clinical trial context.
}

\usage{
Bai.Hu.Shen.Urn(k, p, ssn, Y0 = NULL, nsim = 2000, alpha = 0.05)
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
Bai, Hu and Shen's urn  can be describe as follows: An urn contains \eqn{K} types of balls initially. Balls of types \eqn{1, 2, \cdots, K} represent treatments  \eqn{1, 2, \cdots, K}. A type \eqn{k} ball is drawn randomly from the urn, and then we assign the patient, who is waiting to be assigned, to the treatment \eqn{k}. After obtaining the response, we may adapt the composition of the urn. A success on treatment \eqn{k} adds a ball of type of $k$ to the urn and a failure on treatment \eqn{k} adds \eqn{\frac{p_i}{(M-p_k)}} ball for each of the other \eqn{K-1} types, where \eqn{M = p_1 + ... + p_K}.
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
Bai Z D, Hu F, Shen L. \emph{An adaptive design for multi-arm clinical trials}[J]. Journal of Multivariate Analysis, 2002, 81(1): 1-18.
}

\examples{
## a simple use
bhs.res = Bai.Hu.Shen.Urn(k = 3,
                          p = c(0.7, 0.8, 0.6),
                          ssn = 500,
                          Y0 = NULL,
                          nsim = 200,
                          alpha = 0.05)

## view the output
bhs.res
\donttest{
  ## view all simulation settings
  bhs.res$name
  bhs.res$parameter
  bhs.res$k

  ## View the simulations results
  bhs.res$propotion
  bhs.res$failRate
  bhs.res$pwCalc
  bhs.res$assignment
  }
}
\keyword{BaiHuShenUrn}
\keyword{Bai.Hu.Shen.Urn}
