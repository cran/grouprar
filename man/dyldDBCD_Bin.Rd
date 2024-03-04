\name{dyldDBCD_Bin}
\alias{dyldDBCD_Bin}
\title{
Hu and Zhang's Doubly Biased Coin Deisgn with delayed Binary Response
}
\description{
Simulating Hu and Zhang's doubly biased coin deisgn with delayed binary response (number of arms \eqn{\ge 2}) in a clinical trial context. (Inference: two-sided hypothesis testing t-test or chi-square test)
}
\usage{
dyldDBCD_Bin(n0 = 20, p, k, ssn, ent.param, rspT.dist,
             rspT.param, theta0 = NULL, target.alloc = "RPW", r = 2,
             nsim = 2000, mRate = NULL, alpha = 0.05)
}
\arguments{
  \item{n0}{A positive integer. \code{n0} represents the initial patient population assogned through restricted randomization for initial parameter estimation.}
  \item{p}{A positive vector of length equals to \code{k}. The values specify the true success rates for the various treatments, and these rates are used to generate data for simulations.}
  \item{k}{A positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k = 2})}
  \item{ssn}{A positive integer. The value specifies the total number of participants involved in each round of the simulation.}
  \item{ent.param}{A positive integer. The value specified the parameter for an expoential distribution which determine the time for each participant enter the trial.}
  \item{rspT.dist}{Distribution Type. Specifies the type of distribution that models the time spent for the availability of patient \eqn{i} under treatment \eqn{k}. Acceptable options for this argument include: \code{"exponential"}, \code{"normal"}, and \code{"uniform"}.}
  \item{rspT.param}{A vector. Specifies the parameters required by the distribution that models the time spent for the availability under each treatment. (eg. If there are 3 treatments groups and each of them follows truncated normal distribution with parameter pair (3, 2), (2, 1), (4, 1), repectively. Then the \code{rspT.param = c(3, 2, 2, 1, 4, 1)})}
  \item{theta0}{A vector of length k. Each value in the vector represents a probability used for adjusting parameter estimates. If the argument is not provided, it defaults to a vector of length k, with all values set to 0.5.}
  \item{target.alloc}{Desired allocation proportion. The option for this argument could be one of \code{"Neyman"}, \code{"RSIHR"}, \code{"RPW"}, \code{"WeisUrn"}. The default is \code{"RPW"}. }
  \item{r}{A positive number. Parameter for Hu and Zhang's doubly biased coin design and usually take values 2-4. The default value is 2.}
  \item{nsim}{a positive integer. The value specifies the total number of simulations, with a default value of 2000.}
  \item{mRate}{a numerical value between 0 and 1, inclusive, representing the missing rate for the responses. This parameter pertains to missing-at-random data. The default value is \code{NULL}, indicating no missing values by default.}
  \item{alpha}{a numerical value between 0 and 1. The value represents the predetermined level of significance that defines the probability threshold for rejecting the null hypothesis, with a default value of 0.05.}

}
\details{Hu and Zhang's Doubly Biased Coin Design with delayed Binary Response employs the following treatment allocation scheme:

(a) Initially, due to limited information about treatment efficacy, the first \code{n0} patients are assigned to K treatments using restricted randomization (as described by Rosenberger and Lachin, 2002).

(b) For \eqn{m \ge n_0}, patient \eqn{(m+1)} is allocated to treatment \eqn{k} with a probability \eqn{p_{m+1, k}}, which depends on the available responses and estimated target allocation via \eqn{g_k}, as proposed by Hu and Zhang (2004).

For a more comprehensive description of the procedure, please refer to the paper titled 'Doubly adaptive biased coin designs with delayed responses' authored by Hu et al. in 2008.}

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
Hu, F., Zhang, L. X., Cheung, S. H., & Chan, W. S. (2008). \emph{Doubly adaptive biased coin designs with delayed responses.} Canadian Journal of Statistics, 36(4), 541-559.
}

\examples{
# a simple use
# Define the arguments
## Arguments for generate the simulated data
### For response simulation
p = c(0.6, 0.8)
k = 2
ssn = 100
### for enter time and response time simulation
ent.param = 0.7
rspT.dist = "exponential"
rspT.param = c(1, 1, 3, 1)

## Arguments for the deisgn
n0 = 20
target.alloc = "RSIHR"

dyldDBCD_Bin(n0 = n0, p = p, k = k, ssn = ssn, ent.param, rspT.dist, rspT.param, theta0 = NULL,
             target.alloc, r = 2, nsim = 150, mRate = NULL, alpha = 0.05)
}

\keyword{dyldDBCD_Bin}
\keyword{Doubly adaptive biased coin designs with delayed responses}
