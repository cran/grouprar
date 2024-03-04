\name{Group.DBCD_Bin}
\alias{Group.DBCD_Bin}
\title{
Group Doubly Biased Coin Design with Binary Response Type
}
\description{
Simulating group doubly biased coin deisgn with binary response (number of arms \eqn{\ge 2}) in a clinical trial context. (Inference: two-sided hypothesis testing t-test or chi-square test)
}
\usage{
Group.DBCD_Bin(n0 = 20, p, k, gsize.param, ssn, theta0 = NULL,
                target.alloc = "RPW", r = 2, nsim = 2000,
                mRate = NULL, alpha = 0.05)
}
\arguments{
  \item{n0}{A positive integer. \code{n0} represents the initial patient population assigned through restricted randomization for initial parameter estimation.}
  \item{p}{A positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.}
  \item{k}{A positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})}
  \item{gsize.param}{A positive integer. It represents the expected number of people enrolling in a specified interval, assuming that the enrollment rate per unit time follows a Poisson distribution.}
  \item{ssn}{A positive integer. The value specifies the total number of participants involved in each round of the simulation.}
  \item{theta0}{A vector of length k. Each value in the vector represents a probability used for adjusting parameter estimates. If the argument is not provided, it defaults to a vector of length k, with all values set to 0.5.}
  \item{target.alloc}{Desired allocation proportion. The option for this argument could be one of \code{"Neyman"}, \code{"RSIHR"}, \code{"RPW"}, \code{"WeisUrn"}. The default is \code{"RPW"}. }
  \item{r}{A positive number. Parameter for Hu and Zhang's doubly biased coin design and usually take values 2-4. The default value is 2.}
  \item{nsim}{A positive integer. The value specifies the total number of simulations, with a default value of 2000.}
   \item{mRate}{A numerical value between 0 and 1, inclusive, representing the missing rate for the responses. This parameter pertains to missing-at-random data. The default value is \code{NULL}, indicating no missing values by default.}
  \item{alpha}{A number between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.}

}
\details{
Hu and Zhang's Doubly Biased Coin Design (DBCD) adjusts the probability of assigning each patient to a specific treatment group in a clinical trial, based on the responses of all previous patients. The Group DBCD is an enhanced version of this approach, offering a more practical perspective. It dynamically updates the allocation probabilities for patients in each group based on the responses of all preceding groups, either when the data available or at fixed time intervals (weekly or biweekly).

The process begins by assigning \code{n0} (perhaps the first few groups of) patients to treatment groups using restricted randomization and collecting their responses. Initial parameter estimates for the response variable are then obtained for each treatment group. Subsequently, based on these parameter estimates, the esmatied desired allocation proportion is calculated. Afterward, the Hu and Zhang's allocation function is applied to determine the probabilities for the next group of patients to be assigned to each treatment group, which force the allocation proportion close to the desired one. This process is repeated sequentially for each group until the desired number of patients has been allocated, as predetermined.

This methodology was introduced by Zhai, Li, Zhang and Hu in their 2023 paper titled 'Group Response-Adaptive Randomization with Delayed and Missing Responses'.
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
Zhai, G., Li, Y., Zhang, L. X. & Hu, F. (2023). \emph{Group Response-Adaptive Randomization with Delayed and Missing Responses.}
}

\examples{
Group.DBCD_Bin(n0 = 20, p = c(0.65, 0.8), k = 2, gsize.param = 5,
               ssn = 300, theta0 = NULL, target.alloc = "RPW",
               r = 2, nsim = 400, mRate = NULL, alpha = 0.05)
}

\keyword{group doubly biased coin deisgn}
\keyword{Group.DBCD_Bin}
\keyword{groupDBCD}
\keyword{group doubly biased coin deisgn with discrete response}
