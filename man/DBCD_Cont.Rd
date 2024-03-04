\name{DBCD_Cont}
\alias{DBCD_Cont}
\title{
Hu and Zhang's Doubly Biased Coin Design with Continuous Response Type
}
\description{
Simulating Hu and Zhang's doubly biased coin deisgn with continuous response (number of arms \eqn{\ge 2}) in a clinical trial context. (Inference: two-sided hypothesis testing t-test or chi-square test)
}
\usage{
DBCD_Cont(n0 = 20, theta, k, ssn, theta0 = NULL, target.alloc = "Neyman",
          r = 2, nsim = 2000, alpha = 0.05)
}
\arguments{
  \item{n0}{A positive integer. \code{n0} represents the initial patient population assogned through restricted randomization for initial parameter estimation.}
  \item{theta}{A numerical vector of length equal to \code{2k}. These values specify the true parameters for each treatment and are used for generating data in simulations. For example, if \code{k=2}, you should provide two pairs of parameter values, each consisting of the mean and variance, like: \code{theta = c(13, 4.0^2, 15, 2.5^2)}.}
  \item{k}{A positive integer. The value specifies the number of treatment groups involved in a clinical trial. (\eqn{k \ge 2})}
  \item{ssn}{A positive integer. The value specifies the total number of participants involved in each round of the simulation.}
  \item{theta0}{A vector of length 2k. Each value in the vector represents a probability used for adjusting parameter estimates. If the argument is not provided, it defaults to a vector of length 2k, with all paramter pair setted to be (0, 1).}
  \item{target.alloc}{Desired allocation proportion. The option for this argument could be one of \code{"Neyman"}, \code{"ZR"}, \code{"DaOptimal"}. The default is \code{"Neyman"}. The details see Zhang L. and Rosenberger. W (2006).}
  \item{r}{A positive number. Parameter for Hu and Zhang's doubly biased coin design and usually take values 2-4. The default value is 2.}
  \item{nsim}{a positive integer. The value specifies the total number of simulations, with a default value of 2000.}
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
\item{failRate}{The average response value for the entire trial.}
\item{pwClac}{The probability of the study to detect a significant difference or effect if it truly exists.}
\item{k}{Number of arms involved in the trial.}
}

\references{
Hu, F., Zhang, L. X., Cheung, S. H., & Chan, W. S. (2008). \emph{Doubly adaptive biased coin designs with delayed responses.} Canadian Journal of Statistics, 36(4), 541-559.

Zhang, L., & Rosenberger, W. F. (2006). \emph{Response‐adaptive randomization for clinical trials with continuous outcomes.} Biometrics, 62(2), 562-569.
}

\seealso{
 See \code{\link{DBCD_Bin}} for simulations of Hu and Zhang's doubly biased coin deisgn with binary response.

 See \code{\link{dyldDBCD_Cont}} for simulations of Hu and Zhang's doubly biased coin deisgn with delayed continuous response.
}
\examples{
# A simple use!
# Define the arguments
## Arguments for generate the simulated data
theta = c(13, 4.0^2, 15, 2.5^2)
k = 2
ssn = 88
### Other arguments
target.alloc = "Neyman"

res = DBCD_Cont(n0 = 20, theta, k, ssn, theta0 = NULL,
                target.alloc = "Neyman", r = 2, nsim = 200, alpha = 0.05)

# View the output (A list of all results)
res
}
\keyword{doubly biased coin deisgn}
\keyword{DBCD_Bin}
\keyword{doubly biased coin deisgn with contiunous response}

