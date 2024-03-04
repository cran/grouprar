\name{Group.dyldDBCD_Cont}
\alias{Group.dyldDBCD_Cont}
\title{
Group Doubly Biased Coin Design with Delayed Continuous Response
}
\description{
Simulating group doubly biased coin deisgn with delayed continuous response (number of arms \eqn{\ge 2}) in a clinical trial context. (Inference: two-sided hypothesis testing t-test or chi-square test)
}
\usage{
Group.dyldDBCD_Cont(n0 = 20, theta, k, ssn,
                    gsize.param, rspT.dist, rspT.param,
                    target.alloc = "Neyman",  r = 2,
                    nsim = 2000, eTime = 7,
                    mRate = NULL, alpha = 0.05)
}
\arguments{
\item{n0}{A positive integer. \code{n0} represents the initial patient population assogned through restricted randomization for initial parameter estimation.}
  \item{theta}{A numerical vector of length equal to \code{2k}. These values specify the true parameters for each treatment and are used for generating data in simulations. For example, if \code{k=2}, you should provide two pairs of parameter values, each consisting of the mean and variance, like: \code{theta = c(13, 4.0^2, 15, 2.5^2)}.}
  \item{k}{A positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k = 2})}
  \item{ssn}{A positive integer. The value specifies the total number of participants involved in each round of the simulation.}
  \item{gsize.param}{A positive integer. It represents the expected number of people enrolling in a specified interval, assuming that the enrollment rate per unit time follows a Poisson distribution.}
  \item{rspT.dist}{Distribution Type. Specifies the type of distribution that models the time spent for the availability of patient \eqn{i} under treatment \eqn{k}. Acceptable options for this argument include: \code{"exponential"}, \code{"normal"}, and \code{"uniform"}.}
  \item{rspT.param}{A vector. Specifies the parameters required by the distribution that models the time spent for the availability under each treatment. (eg. If there are 3 treatments groups and each of them follows truncated normal distribution with parameter pair (3, 2), (2, 1), (4, 1), repectively. Then the \code{rspT.param = c(3, 2, 2, 1, 4, 1)})}
  \item{target.alloc}{Desired allocation proportion. The option for this argument could be one of \code{"Neyman"}, \code{"ZR"}, \code{"DaOptimal"}. The default is \code{"Neyman"}. The details see Zhang L. and Rosenberger. W (2006).}
  \item{r}{A positive number. Parameter for Hu and Zhang's doubly biased coin design and usually take values 2-4. The default value is 2.}
  \item{nsim}{a positive integer. The value specifies the total number of simulations, with a default value of 2000.}
  \item{eTime}{A positive number. The interval time between enrollment of participants in each group. The default is 7.}
  \item{mRate}{a numerical value between 0 and 1, inclusive, representing the missing rate for the responses. This parameter pertains to missing-at-random data. The default value is \code{NULL}, indicating no missing values by default.}
  \item{alpha}{a numerical value between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.}
}


\details{
Hu and Zhang's Doubly Biased Coin Design (DBCD) adjusts the probability of assigning each patient to a specific treatment group in a clinical trial, based on the responses of all previous patients. The Group DBCD is an enhanced version of this approach, offering a more practical perspective. It dynamically updates the allocation probabilities for patients in each group based on the avaiable responses of all preceding groups, either when the data available or at fixed time intervals (weekly or biweekly). Here, the function focuses on implementing the group doubly biased coin design, tailored to delayed continuous responses.

The process begins by assigning \code{n0} (perhaps the first few groups of) patients to treatment groups using restricted randomization and collecting their responses. Initial parameter estimates for the response variable are then obtained for each treatment group. Subsequently, based on these parameter estimates, the esmatied desired allocation proportion is calculated. Afterward, the Hu and Zhang's allocation function is applied to determine the probabilities for the next group of patients to be assigned to each treatment group, which force the allocation proportion close to the desired one. This process is repeated sequentially for each group until the desired number of patients has been allocated, as predetermined.

This methodology was introduced by Zhai, Li, Zhang and Hu in their 2023 paper titled 'Group Response-Adaptive Randomization with Delayed and Missing Responses'.
}
\value{
\item{name}{The name of procedure.}
\item{parameter}{The true parameters used to do the simulations.}
\item{assignment}{The randomization sequence.}
\item{propotion}{Average allocation porpotion for each of treatment groups.}
\item{failRate}{The average response value for the entire trial.}
\item{pwClac}{The probability of the study to detect a significant difference or effect if it truly exists.}
\item{k}{Number of arms involved in the trial.}
}
\references{
Zhai, G., Li, Y., Zhang, L. X. & Hu, F. (2023). \emph{Group Response-Adaptive Randomization with Delayed and Missing Responses.}
}

\examples{
# a simple use
# Define the arguments
## Arguments for generate the simulated data
### For response simulation
theta = c(13, 4.0^2, 15, 2.5^2)
k = 2
ssn = 120

### for enter time and response time simulation
eTime = 7
gsize.param = 5
rspT.param = rep(10, 2)
rspT.dist = "exponential"

## Arguments for the deisgn
n0 = 10
target.alloc = "Neyman"

res = Group.dyldDBCD_Cont(n0, theta, k, ssn, gsize.param,
                          rspT.dist, rspT.param, target.alloc = "Neyman",
                          r = 2, nsim = 200, eTime = 7,  mRate = 0.2, alpha = 0.05)

# View the output (A list of all results)
res
}

\keyword{group doubly biased coin deisgn}
\keyword{Group.dyldDBCD_Cont}
\keyword{groupDBCD}
\keyword{group doubly biased coin deisgn with delayed continuous response}

