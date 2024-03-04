\name{DBCD_Bin}
\alias{DBCD_Bin}
\title{
Hu and Zhang's Doubly Biased Coin Design with Binary Response Type
}
\description{
Simulating Hu and Zhang's doubly biased coin deisgn with binary response (number of arms \eqn{\ge 2}) in a clinical trial context. (Inference: two-sided hypothesis testing t-test or chi-square test)
}
\usage{
DBCD_Bin(n0 = 20, p, k, ssn, theta0 = NULL, target.alloc = "RPW",
         r = 2, nsim = 2000, mRate = NULL, alpha = 0.05)
}
\arguments{
  \item{n0}{A positive integer. \code{n0} represents the initial patient population assogned through restricted randomization for initial parameter estimation.}
  \item{p}{A positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.}
  \item{k}{A positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})}
  \item{ssn}{A positive integer. The value specifies the total number of participants involved in each round of the simulation.}
  \item{theta0}{A vector of length k. Each value in the vector represents a probability used for adjusting parameter estimates. If the argument is not provided, it defaults to a vector of length k, with all values set to 0.5.}
  \item{target.alloc}{Desired allocation proportion. The option for this argument could be one of \code{"Neyman"}, \code{"RSIHR"}, \code{"RPW"}, \code{"WeisUrn"}. The default is \code{"RPW"}. }
  \item{r}{A positive number. Parameter for Hu and Zhang's doubly biased coin design and usually take values 2-4. The default value is 2.}
  \item{nsim}{a positive integer. The value specifies the total number of simulations, with a default value of 2000.}
  \item{mRate}{a numerical value between 0 and 1, inclusive, representing the missing rate for the responses. This parameter pertains to missing-at-random data. The default value is \code{NULL}, indicating no missing values by default.}
  \item{alpha}{A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.}

}
\details{
The objective of Hu and Zhang's doubly biased coin design is to allocate patients sequentially while closely approximating the desired allocation proportion, which is a function of certain unknown parameters related to the response variable under each treatment.

The process begins by assigning \code{n0} patients to treatment groups using restricted randomization and collecting their responses. Initial parameter estimates for the response variable are then obtained for each treatment group. Subsequently, based on these parameter estimates, the desired allocation proportion is calculated. Afterward, the Hu and Zhang's allocation function is applied to determine the probabilities for the next patient to be assigned to each treatment group, which force the allocation proportion close to the desired one. This process is repeated sequentially for each patient until the desired number of patients has been allocated, as predetermined.

This methodology was introduced by Hu and Zhang in their 2004 paper titled 'Asymptotic Properties of Doubly Adaptive Biased Coin Designs for Multitreatment Clinical Trials.
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
Hu, F., & Zhang, L. X. (2004). \emph{Asymptotic properties of doubly adaptive biased coin designs for multi-treatment clinical trials.} The Annals of Statistics, 32(1), 268-301.
}

\examples{
DBCD_Bin(n0 = 20, p = c(0.7, 0.8), k = 2, ssn = 300, theta0 = NULL,
         target.alloc = "RPW", r = 2, nsim = 50, mRate = NULL, alpha = 0.05)
}

\keyword{doubly biased coin deisgn}
\keyword{DBCD_Bin}
