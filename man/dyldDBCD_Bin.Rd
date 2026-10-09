\name{dyldDBCD_Bin}
\alias{dyldDBCD_Bin}
\title{
Hu and Zhang's Doubly Biased Coin Design with Delayed Binary Response
}
\description{
Simulating Hu and Zhang's doubly biased coin design with delayed binary response (number of arms \eqn{\ge 2}) in a clinical trial context. (Inference: two-sided t-test for two arms or chi-square test for more than two arms.)
}
\usage{
dyldDBCD_Bin(n0 = 20, p, k, ssn, ent.param, rspT.dist, rspT.param,
             theta0 = NULL, target.alloc = "RPW", r = 2, nsim = 2000,
             mRate = NULL, alpha = 0.05, allocation = "DBCD",
             erade.alpha = 0.5, lower.bound = 0, test.fun = NULL,
             typeI = FALSE, seed = NULL)
}
\arguments{
  \item{n0}{A positive integer and a multiple of \code{k}. The number of initial patients assigned by restricted randomization for initial parameter estimation.}
  \item{p}{A vector of length \code{k} with values between 0 and 1. The true success rates of the treatments, used to generate data for the simulations.}
  \item{k}{A positive integer. The number of treatment groups in the trial (\eqn{k \ge 2}).}
  \item{ssn}{A positive integer. The total number of participants in each simulated trial.}
  \item{ent.param}{A positive number. The mean time between consecutive patient arrivals. Inter-arrival times are drawn from an exponential distribution with this mean.}
  \item{rspT.dist}{The distribution of the time from enrollment until the response is observed. One of \code{"exponential"}, \code{"normal"} or \code{"uniform"}.}
  \item{rspT.param}{A numeric vector of parameters for the response-time distributions, one distribution per treatment and response, ordered as (treatment 1 failure, treatment 1 success, treatment 2 failure, treatment 2 success, ...). For \code{"exponential"} give the \eqn{2k} means. For \code{"normal"} give \eqn{2k} (mean, sd) pairs and for \code{"uniform"} give \eqn{2k} (min, max) pairs, that is \eqn{4k} values in total. Normal times below 0 are set to 0. For example, with \code{k = 2} and exponential times, \code{rspT.param = c(3, 2, 4, 1)} gives mean times 3 and 2 for failures and successes on treatment 1, and 4 and 1 on treatment 2.}
  \item{theta0}{A vector of length \code{k} used to smooth the success-rate estimates, \eqn{\hat p_k = (S_k + \theta_{0k})/(N_k + 1)}{p_k = (S_k + theta0_k)/(N_k + 1)}, where \eqn{S_k} and \eqn{N_k} are the number of observed successes and observed responses on treatment \eqn{k}. If \code{NULL} (default), all values are 0.5.}
  \item{target.alloc}{Desired allocation proportion. One of \code{"Neyman"}, \code{"RSIHR"}, \code{"RPW"}, \code{"WeisUrn"}, \code{"OptimalNeyman"} or \code{"OptimalRSIHR"}. The default is \code{"RPW"}. See Details.}
  \item{r}{A non-negative number. The tuning parameter of Hu and Zhang's allocation function, used when \code{allocation = "DBCD"}. Values between 2 and 4 are common. The default value is 2.}
  \item{nsim}{A positive integer. The number of simulated trials, with a default value of 2000.}
  \item{mRate}{A number between 0 and 1 giving the probability that a response is missing. Missingness is simulated completely at random (MCAR), and missing responses are never observed and are excluded from estimation and testing. The default \code{NULL} means no missing responses. As in Zhai et al. (2024), missing responses are excluded from the estimates, the failure rate and the test, while the allocation proportions (both those used by the allocation function and the reported ones) count all enrolled patients.}
  \item{alpha}{A number between 0 and 1. The significance level of the two-sided test, with a default value of 0.05.}
  \item{allocation}{The allocation function. \code{"DBCD"} (default) uses the allocation function of Hu and Zhang (2004) with parameter \code{r}. \code{"ERADE"} uses the efficient randomized adaptive design of Hu, Zhang and He (2009), in its multi-arm version (Alkhnefr, Hu and Zhai, 2025), with parameter \code{erade.alpha}. See Details.}
  \item{erade.alpha}{A number between 0 and 1. The degree of randomization of ERADE, used when \code{allocation = "ERADE"}. Smaller values push the allocation more strongly towards the target, and values of 1/2 or 2/3 are recommended. The default is 0.5.}
  \item{lower.bound}{A number between 0 and \eqn{1/k}. The smallest allowed target proportion of each arm for the targets \code{"OptimalNeyman"} and \code{"OptimalRSIHR"}. The default is 0.}
  \item{test.fun}{An optional function \code{function(outcome, assignment)} that returns the p-value of a test for one simulated trial, where \code{outcome} holds the observed responses and \code{assignment} the treatment labels (integers 1 to \code{k}); missing responses are removed first. The null hypothesis is rejected when the p-value is at most \code{alpha}. The default \code{NULL} uses the built-in two-sided tests: the t-test for two arms and the Wald chi-squared test of equal means for more than two arms (for binary responses, if some arms have no variability, the variances are computed from \eqn{(S_k + 1)/(n_k + 2)}{(S_k + 1)/(n_k + 2)} as in Agresti and Caffo, 2000).}
  \item{typeI}{Logical. If \code{TRUE}, the simulation is repeated under the null hypothesis, with every arm's success rate set to the average of \code{p}, and the rejection rate is reported as \code{type I error}. The default is \code{FALSE}.}
  \item{seed}{An optional integer passed to \code{\link{set.seed}} before simulating, so that the results can be reproduced. The default \code{NULL} leaves the random number generator unchanged.}
}
\details{Hu and Zhang's doubly biased coin design with delayed binary response uses the following allocation scheme.

(a) Initially, due to limited information about treatment efficacy, the first \code{n0} patients are assigned to the K treatments using restricted randomization (as described by Rosenberger and Lachin, 2002).

(b) For \eqn{m \ge n_0}, patient \eqn{(m+1)} is allocated to treatment \eqn{k} with probability \eqn{p_{m+1, k}}, which depends on the responses available at that time and on the estimated target allocation through the allocation function \eqn{g_k} proposed by Hu and Zhang (2004).

For a more comprehensive description of the procedure, please refer to the paper 'Doubly adaptive biased coin designs with delayed responses' by Hu et al. (2008).

Target allocations. With \eqn{q_k = 1 - p_k}, \code{"Neyman"} and \code{"RSIHR"} are proportional to \eqn{\sqrt{p_k q_k}}{sqrt(p_k q_k)} and \eqn{\sqrt{p_k}}{sqrt(p_k)}. For two arms these are the Neyman allocation and the optimal allocation of Rosenberger et al. (2001), and for more arms they are simple generalizations. \code{"RPW"} and \code{"WeisUrn"} are proportional to \eqn{1/q_k}, the limiting allocation of the randomized play-the-winner rule and of Wei's urn. \code{"OptimalNeyman"} and \code{"OptimalRSIHR"} are the k-arm optimal allocations of Tymofyeyev, Rosenberger and Hu (2007). They minimize the total sample size and the expected number of failures, respectively, for a fixed noncentrality parameter of the chi-squared test of equal success rates, with every proportion at least \code{lower.bound}. For two arms and \code{lower.bound = 0} they equal \code{"Neyman"} and \code{"RSIHR"}. For three or more arms the optimum without a lower bound can assign no patients to the middle arms, so a positive \code{lower.bound} such as 0.1 is recommended.

Allocation functions. With \code{allocation = "ERADE"} and \eqn{\alpha} = \code{erade.alpha}, arm \eqn{k} receives the next patient (in group designs, each patient of the next group) with probability \eqn{\alpha\hat\rho_k}{alpha * rho_k} if it is over-allocated (\eqn{N_k/m > \hat\rho_k}{N_k/m > rho_k}), \eqn{\hat\rho_k}{rho_k} if it is on target, and \eqn{\hat\rho_k + (1-\alpha)\sum_{j \in S}\hat\rho_j/|T|}{rho_k + (1 - alpha) sum_{j in S} rho_j / |T|} if it is under-allocated, where \eqn{S} and \eqn{T} are the sets of over- and under-allocated arms, \eqn{N_k/m} is the current allocation proportion and \eqn{\hat\rho}{rho} is the estimated target (Alkhnefr, Hu and Zhai, 2025, eq. (1)). For two arms this is the ERADE of Hu, Zhang and He (2009). ERADE attains the lower bound of the asymptotic variance of the allocation proportions, so its allocation is less variable than that of the DBCD.
}

\value{
An object of class \code{"grouprar"}, a list that is printed as a short summary (see \code{\link{print.grouprar}}), with the following elements.
\item{method}{The name of the procedure.}
\item{sample size}{The total sample size (\code{Total}) and the expected number of patients with observed responses (\code{Effective Size}).}
\item{parameter}{The true success rates used in the simulations, named \code{pA}, \code{pB}, ...}
\item{propotion}{The mean allocation proportion of each arm over the simulations, named \code{treatment A}, \code{treatment B}, ...}
\item{sd of propotion}{The standard deviation of the allocation proportion of each arm over the simulations.}
\item{failure rate}{The mean failure rate over the simulations.}
\item{sd of failure rate}{The standard deviation of the failure rate (or mean response) over the simulations.}
\item{power}{The proportion of simulated trials that reject the null hypothesis of equal success rates. Simulations in which the test cannot be computed are dropped.}
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

Rosenberger, W. F., Stallard, N., Ivanova, A., Harper, C. N. and Ricks, M. L. (2001). \emph{Optimal adaptive designs for binary response trials.} Biometrics, 57(3), 909-913. \doi{10.1111/j.0006-341X.2001.00909.x}

Tymofyeyev, Y., Rosenberger, W. F. and Hu, F. (2007). \emph{Implementing optimal allocation in sequential binary response experiments.} Journal of the American Statistical Association, 102(477), 224-234. \doi{10.1198/016214506000000906}
}

\examples{
# a simple use
## Arguments for generating the simulated data
### For response simulation
p = c(0.6, 0.8)
k = 2
ssn = 100
### for entry time and response time simulation
ent.param = 0.7
rspT.dist = "exponential"
rspT.param = c(1, 1, 3, 1)

## Arguments for the design
n0 = 20
target.alloc = "RSIHR"

res <- dyldDBCD_Bin(n0 = n0, p = p, k = k, ssn = ssn, ent.param = ent.param,
                    rspT.dist = rspT.dist, rspT.param = rspT.param, theta0 = NULL,
                    target.alloc = target.alloc, r = 2, nsim = 100,
                    mRate = NULL, alpha = 0.05)
res
}

\keyword{design}
\concept{doubly adaptive biased coin design}
\concept{doubly adaptive biased coin designs with delayed responses}
\concept{response-adaptive randomization}
