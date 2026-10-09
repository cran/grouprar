\name{dyldDBCD_Cont}
\alias{dyldDBCD_Cont}
\title{
Hu and Zhang's Doubly Biased Coin Design with Delayed Continuous Response
}
\description{
Simulating Hu and Zhang's doubly biased coin design with delayed continuous response (number of arms \eqn{\ge 2}) in a clinical trial context. (Inference: two-sided t-test for two arms or chi-square test for more than two arms.)
}
\usage{
dyldDBCD_Cont(n0 = 20, theta, k, ssn, ent.param, rspT.dist, rspT.param,
              target.alloc = "Neyman", r = 2, nsim = 2000, mRate = NULL,
              alpha = 0.05, allocation = "DBCD", erade.alpha = 0.5,
              lower.bound = 0, test.fun = NULL, typeI = FALSE, seed = NULL)
}
\arguments{
  \item{n0}{A positive integer and a multiple of \code{k}. The number of initial patients assigned by restricted randomization for initial parameter estimation.}
  \item{theta}{A numerical vector of length \code{2k} giving the true mean and variance of each treatment, used to generate data for the simulations. For example, with \code{k = 2} use \code{theta = c(13, 4.0^2, 15, 2.5^2)}.}
  \item{k}{A positive integer. The number of treatment groups in the trial (\eqn{k \ge 2}).}
  \item{ssn}{A positive integer. The total number of participants in each simulated trial.}
  \item{ent.param}{A positive number. The mean time between consecutive patient arrivals. Inter-arrival times are drawn from an exponential distribution with this mean.}
  \item{rspT.dist}{The distribution of the time from enrollment until the response is observed. One of \code{"exponential"}, \code{"normal"} or \code{"uniform"}.}
  \item{rspT.param}{A numeric vector of parameters for the response-time distribution of each treatment. For \code{"exponential"} give the \eqn{k} means. For \code{"normal"} give \eqn{k} (mean, sd) pairs and for \code{"uniform"} give \eqn{k} (min, max) pairs, that is \eqn{2k} values in total. Normal times below 0 are set to 0. For example, with 3 treatments and normal times with (mean, sd) pairs (3, 2), (2, 1) and (4, 1), use \code{rspT.param = c(3, 2, 2, 1, 4, 1)}.}
  \item{target.alloc}{Desired allocation proportion. One of \code{"Neyman"}, \code{"ZR"}, \code{"OptimalNeyman"} or \code{"DaOptimal"}. The default is \code{"Neyman"}. See Details.}
  \item{r}{A non-negative number. The tuning parameter of Hu and Zhang's allocation function, used when \code{allocation = "DBCD"}. Values between 2 and 4 are common. The default value is 2.}
  \item{nsim}{A positive integer. The number of simulated trials, with a default value of 2000.}
  \item{mRate}{A number between 0 and 1 giving the probability that a response is missing. Missingness is simulated completely at random (MCAR), and missing responses are never observed and are excluded from estimation and testing. The default \code{NULL} means no missing responses. As in Zhai et al. (2024), missing responses are excluded from the estimates, the failure rate and the test, while the allocation proportions (both those used by the allocation function and the reported ones) count all enrolled patients.}
  \item{alpha}{A number between 0 and 1. The significance level of the two-sided test, with a default value of 0.05.}
  \item{allocation}{The allocation function. \code{"DBCD"} (default) uses the allocation function of Hu and Zhang (2004) with parameter \code{r}. \code{"ERADE"} uses the efficient randomized adaptive design of Hu, Zhang and He (2009), in its multi-arm version (Alkhnefr, Hu and Zhai, 2025), with parameter \code{erade.alpha}. See Details.}
  \item{erade.alpha}{A number between 0 and 1. The degree of randomization of ERADE, used when \code{allocation = "ERADE"}. Smaller values push the allocation more strongly towards the target, and values of 1/2 or 2/3 are recommended. The default is 0.5.}
  \item{lower.bound}{A number between 0 and \eqn{1/k}. The smallest allowed target proportion of each arm for the targets \code{"ZR"} and \code{"OptimalNeyman"}. The default is 0.}
  \item{test.fun}{An optional function \code{function(outcome, assignment)} that returns the p-value of a test for one simulated trial, where \code{outcome} holds the observed responses and \code{assignment} the treatment labels (integers 1 to \code{k}); missing responses are removed first. The null hypothesis is rejected when the p-value is at most \code{alpha}. The default \code{NULL} uses the built-in two-sided tests: the t-test for two arms and the Wald chi-squared test of equal means for more than two arms (for binary responses, if some arms have no variability, the variances are computed from \eqn{(S_k + 1)/(n_k + 2)}{(S_k + 1)/(n_k + 2)} as in Agresti and Caffo, 2000).}
  \item{typeI}{Logical. If \code{TRUE}, the simulation is repeated under the null hypothesis, with the mean of every arm set to the average of the means in \code{theta} and the variance of each arm kept, and the rejection rate is reported as \code{type I error}. The default is \code{FALSE}.}
  \item{seed}{An optional integer passed to \code{\link{set.seed}} before simulating, so that the results can be reproduced. The default \code{NULL} leaves the random number generator unchanged.}
}
\details{Hu and Zhang's doubly biased coin design with delayed continuous response uses the following allocation scheme.

(a) Initially, due to limited information about treatment efficacy, the first \code{n0} patients are assigned to the K treatments using restricted randomization (as described by Rosenberger and Lachin, 2002).

(b) For \eqn{m \ge n_0}, patient \eqn{(m+1)} is allocated to treatment \eqn{k} with probability \eqn{p_{m+1, k}}, which depends on the responses available at that time and on the estimated target allocation through the allocation function \eqn{g_k} proposed by Hu and Zhang (2004).

For a more comprehensive description of the procedure, please refer to the paper 'Doubly adaptive biased coin designs with delayed responses' by Hu et al. (2008).

Target allocations. \code{"Neyman"} is proportional to \eqn{\sigma_k}{sigma_k}. \code{"ZR"} is the allocation of Zhang and Rosenberger (2006), which minimizes the total expected response (smaller responses are better, and all means must be positive). For two arms it is their rule (7): \eqn{\rho_1 = \sigma_1\sqrt{\mu_2}/(\sigma_1\sqrt{\mu_2} + \sigma_2\sqrt{\mu_1})}{rho_1 = sigma_1 sqrt(mu_2)/(sigma_1 sqrt(mu_2) + sigma_2 sqrt(mu_1))} when this assigns more patients to the arm with the smaller mean, and 1/2 otherwise. For three or more arms it minimizes the total expected response for a fixed noncentrality parameter of the chi-squared test of equal means, with every proportion at least \code{lower.bound}, the continuous analogue of Tymofyeyev, Rosenberger and Hu (2007). \code{"OptimalNeyman"} minimizes the total sample size under the same constraints and equals \code{"Neyman"} for two arms with \code{lower.bound = 0}. For three or more arms a positive \code{lower.bound} such as 0.1 is recommended with \code{"ZR"} and \code{"OptimalNeyman"}. \code{"DaOptimal"} is proportional to \eqn{\sigma_k^{4/3}}{sigma_k^(4/3)}. While the target cannot be estimated (for example with fewer than two responses in an arm) equal allocation is used.

Allocation functions. With \code{allocation = "ERADE"} and \eqn{\alpha} = \code{erade.alpha}, arm \eqn{k} receives the next patient (in group designs, each patient of the next group) with probability \eqn{\alpha\hat\rho_k}{alpha * rho_k} if it is over-allocated (\eqn{N_k/m > \hat\rho_k}{N_k/m > rho_k}), \eqn{\hat\rho_k}{rho_k} if it is on target, and \eqn{\hat\rho_k + (1-\alpha)\sum_{j \in S}\hat\rho_j/|T|}{rho_k + (1 - alpha) sum_{j in S} rho_j / |T|} if it is under-allocated, where \eqn{S} and \eqn{T} are the sets of over- and under-allocated arms, \eqn{N_k/m} is the current allocation proportion and \eqn{\hat\rho}{rho} is the estimated target (Alkhnefr, Hu and Zhai, 2025, eq. (1)). For two arms this is the ERADE of Hu, Zhang and He (2009). ERADE attains the lower bound of the asymptotic variance of the allocation proportions, so its allocation is less variable than that of the DBCD.
}

\value{
An object of class \code{"grouprar"}, a list that is printed as a short summary (see \code{\link{print.grouprar}}), with the following elements.
\item{method}{The name of the procedure.}
\item{sample size}{The total sample size (\code{Total}) and the expected number of patients with observed responses (\code{Effective Size}).}
\item{parameter}{The true means and variances used in the simulations, named \code{muA}, \code{sigma2A}, \code{muB}, ...}
\item{propotion}{The mean allocation proportion of each arm over the simulations, named \code{treatment A}, \code{treatment B}, ...}
\item{sd of propotion}{The standard deviation of the allocation proportion of each arm over the simulations.}
\item{failure rate}{The mean response over the simulations (for continuous responses this element holds the mean response, not a failure rate).}
\item{sd of failure rate}{The standard deviation of the failure rate (or mean response) over the simulations.}
\item{power}{The proportion of simulated trials that reject the null hypothesis of equal means. Simulations in which the test cannot be computed are dropped.}
\item{data: failureRate}{The failure rate (or mean response) of each simulated trial.}
\item{data: test}{The test decision of each simulated trial (1 = reject).}
\item{data: assignment}{The treatment assignments of the last simulated trial.}
\item{data: propotion}{A data frame with the allocation proportions of each simulated trial.}
\item{data: allocation}{An \code{nsim} by \code{ssn} matrix with the treatment assignments of every simulated trial.}
\item{type I error}{Only if \code{typeI = TRUE}. The rejection rate under the null hypothesis.}
\item{duration}{The mean time from the first enrollment to the last observed response. Trials without any observed response have duration NA and are left out of the mean.}
\item{sd of duration}{The standard deviation of the duration over the simulations.}
\item{enrollment duration}{The mean time from the first to the last enrollment.}
\item{data: duration}{The duration of each simulated trial.}
\item{data: enrollment}{The enrollment duration of each simulated trial.}
}

\references{
Alkhnefr, N., Hu, F. and Zhai, G. (2025). \emph{Efficient randomized adaptive designs for multi-arm clinical trials.} Statistical Methods in Medical Research, 34(9), 1886-1898. \doi{10.1177/09622802251362644}

Hu, F. and Zhang, L. X. (2004). \emph{Asymptotic properties of doubly adaptive biased coin designs for multitreatment clinical trials.} The Annals of Statistics, 32(1), 268-301. \doi{10.1214/aos/1079120137}

Hu, F., Zhang, L. X. and He, X. (2009). \emph{Efficient randomized-adaptive designs.} The Annals of Statistics, 37(5A), 2543-2560. \doi{10.1214/08-AOS655}

Hu, F., Zhang, L. X., Cheung, S. H. and Chan, W. S. (2008). \emph{Doubly adaptive biased coin designs with delayed responses.} Canadian Journal of Statistics, 36(4), 541-559. \doi{10.1002/cjs.5550360404}

Rosenberger, W. F. and Lachin, J. M. (2002). \emph{Randomization in Clinical Trials: Theory and Practice.} John Wiley & Sons.

Tymofyeyev, Y., Rosenberger, W. F. and Hu, F. (2007). \emph{Implementing optimal allocation in sequential binary response experiments.} Journal of the American Statistical Association, 102(477), 224-234. \doi{10.1198/016214506000000906}

Zhang, L. and Rosenberger, W. F. (2006). \emph{Response-adaptive randomization for clinical trials with continuous outcomes.} Biometrics, 62(2), 562-569. \doi{10.1111/j.1541-0420.2005.00496.x}
}

\examples{
# a simple use
## Arguments for generating the simulated data
### For response simulation
theta = c(13, 4.0^2, 15, 2.5^2)
k = 2
ssn = 88

### for entry time and response time simulation
ent.param = 5
rspT.param = rep(10, 2)
rspT.dist = "exponential"

## Arguments for the design
target.alloc = "Neyman"

res <- dyldDBCD_Cont(n0 = 10, theta = theta, k = k, ssn = ssn, ent.param = ent.param,
                     rspT.dist = rspT.dist, rspT.param = rspT.param,
                     target.alloc = target.alloc, r = 2, nsim = 100,
                     mRate = 0.2, alpha = 0.05)

# View the output (a list of all results)
res
}

\keyword{design}
\concept{doubly adaptive biased coin design}
\concept{doubly adaptive biased coin designs with delayed responses}
\concept{response-adaptive randomization}
