\name{Group.DBCD_Cont}
\alias{Group.DBCD_Cont}
\title{
Group Doubly Biased Coin Design with Continuous Response Type
}
\description{
Simulating the group doubly biased coin design with continuous response (number of arms \eqn{\ge 2}) in a clinical trial context. (Inference: two-sided t-test for two arms or chi-square test for more than two arms.)
}
\usage{
Group.DBCD_Cont(n0 = 20, theta, k, gsize.param, ssn, theta0 = NULL,
                target.alloc = "Neyman", r = 2, nsim = 2000, mRate = NULL,
                alpha = 0.05, allocation = "DBCD", erade.alpha = 0.5,
                lower.bound = 0, monitor = NULL, test.fun = NULL,
                typeI = FALSE, seed = NULL)
}
\arguments{
  \item{n0}{A positive integer and a multiple of \code{k}. Whole groups are assigned by restricted randomization until at least \code{n0} patients are enrolled, for initial parameter estimation.}
  \item{theta}{A numerical vector of length \code{2k} giving the true mean and variance of each treatment, used to generate data for the simulations. For example, with \code{k = 2} use \code{theta = c(13, 4.0^2, 15, 2.5^2)}.}
  \item{k}{A positive integer. The number of treatment groups in the trial (\eqn{k \ge 2}).}
  \item{gsize.param}{A positive number. Group sizes are drawn from a zero-truncated Poisson distribution with rate \code{gsize.param}, so every group has at least one patient.}
  \item{ssn}{A positive integer. The total number of participants in each simulated trial.}
  \item{theta0}{Currently unused. Kept for backward compatibility.}
  \item{target.alloc}{Desired allocation proportion. One of \code{"Neyman"}, \code{"ZR"}, \code{"OptimalNeyman"} or \code{"DaOptimal"}. The default is \code{"Neyman"}. See Details.}
  \item{r}{A non-negative number. The tuning parameter of Hu and Zhang's allocation function, used when \code{allocation = "DBCD"}. Values between 2 and 4 are common. The default value is 2.}
  \item{nsim}{A positive integer. The number of simulated trials, with a default value of 2000.}
  \item{mRate}{A number between 0 and 1 giving the probability that a response is missing. Missingness is simulated completely at random (MCAR), and missing responses are excluded from estimation and testing. The default \code{NULL} means no missing responses. As in Zhai et al. (2024), missing responses are excluded from the estimates, the failure rate and the test, while the allocation proportions (both those used by the allocation function and the reported ones) count all enrolled patients.}
  \item{alpha}{A number between 0 and 1. The significance level of the two-sided test, with a default value of 0.05.}
  \item{allocation}{The allocation function. \code{"DBCD"} (default) uses the allocation function of Hu and Zhang (2004) with parameter \code{r}. \code{"ERADE"} uses the efficient randomized adaptive design of Hu, Zhang and He (2009), in its multi-arm version (Alkhnefr, Hu and Zhai, 2025), with parameter \code{erade.alpha}. See Details.}
  \item{erade.alpha}{A number between 0 and 1. The degree of randomization of ERADE, used when \code{allocation = "ERADE"}. Smaller values push the allocation more strongly towards the target, and values of 1/2 or 2/3 are recommended. The default is 0.5.}
  \item{lower.bound}{A number between 0 and \eqn{1/k}. The smallest allowed target proportion of each arm for the targets \code{"ZR"} and \code{"OptimalNeyman"}. The default is 0.}
  \item{monitor}{An optional object created by \code{\link{sqMonitor}} for group sequential monitoring of a two-arm trial. For each information time \code{t}, a look is taken at the end of the first group in which at least \code{ceiling(t * ssn)} patients are enrolled, and the trial stops for efficacy as soon as the Wald statistic crosses the boundary of \code{\link{sqBoundary}}. The first look must come after the \code{n0} initial patients. Cannot be combined with \code{test.fun}. The default \code{NULL} gives a fixed sample design. See Details.}
  \item{test.fun}{An optional function \code{function(outcome, assignment)} that returns the p-value of a test for one simulated trial, where \code{outcome} holds the observed responses and \code{assignment} the treatment labels (integers 1 to \code{k}); missing responses are removed first. The null hypothesis is rejected when the p-value is at most \code{alpha}. The default \code{NULL} uses the built-in two-sided tests: the t-test for two arms and the Wald chi-squared test of equal means for more than two arms (for binary responses, if some arms have no variability, the variances are computed from \eqn{(S_k + 1)/(n_k + 2)}{(S_k + 1)/(n_k + 2)} as in Agresti and Caffo, 2000). Cannot be combined with \code{monitor}.}
  \item{typeI}{Logical. If \code{TRUE}, the simulation is repeated under the null hypothesis, with the mean of every arm set to the average of the means in \code{theta} and the variance of each arm kept, and the rejection rate is reported as \code{type I error}. The default is \code{FALSE}.}
  \item{seed}{An optional integer passed to \code{\link{set.seed}} before simulating, so that the results can be reproduced. The default \code{NULL} leaves the random number generator unchanged.}
}
\details{
Hu and Zhang's doubly biased coin design (DBCD) adjusts the probability of assigning each patient to a specific treatment group in a clinical trial, based on the responses of all previous patients. The group DBCD is a more practical version of this approach. It updates the allocation probabilities for the patients in each group based on the responses of all preceding groups, either when the data become available or at fixed time intervals (for example weekly or biweekly).

The process begins by assigning the first \code{n0} patients (possibly the first few groups) to treatment groups using restricted randomization and collecting their responses. Initial parameter estimates for the response variable are then obtained for each treatment group. Based on these estimates, the estimated desired allocation proportion is calculated. Then Hu and Zhang's allocation function is applied to determine the probabilities for the next group of patients to be assigned to each treatment group, which drive the allocation proportion towards the desired one. This process is repeated sequentially for each group until the predetermined number of patients has been allocated.

This methodology was introduced by Zhai, Li, Zhang and Hu (2024) in their paper 'Group response-adaptive randomization with delayed and missing responses'.

Target allocations. \code{"Neyman"} is proportional to \eqn{\sigma_k}{sigma_k}. \code{"ZR"} is the allocation of Zhang and Rosenberger (2006), which minimizes the total expected response (smaller responses are better, and all means must be positive). For two arms it is their rule (7): \eqn{\rho_1 = \sigma_1\sqrt{\mu_2}/(\sigma_1\sqrt{\mu_2} + \sigma_2\sqrt{\mu_1})}{rho_1 = sigma_1 sqrt(mu_2)/(sigma_1 sqrt(mu_2) + sigma_2 sqrt(mu_1))} when this assigns more patients to the arm with the smaller mean, and 1/2 otherwise. For three or more arms it minimizes the total expected response for a fixed noncentrality parameter of the chi-squared test of equal means, with every proportion at least \code{lower.bound}, the continuous analogue of Tymofyeyev, Rosenberger and Hu (2007). \code{"OptimalNeyman"} minimizes the total sample size under the same constraints and equals \code{"Neyman"} for two arms with \code{lower.bound = 0}. For three or more arms a positive \code{lower.bound} such as 0.1 is recommended with \code{"ZR"} and \code{"OptimalNeyman"}. \code{"DaOptimal"} is proportional to \eqn{\sigma_k^{4/3}}{sigma_k^(4/3)}. While the target cannot be estimated (for example with fewer than two responses in an arm) equal allocation is used.

Allocation functions. With \code{allocation = "ERADE"} and \eqn{\alpha} = \code{erade.alpha}, arm \eqn{k} receives the next patient (in group designs, each patient of the next group) with probability \eqn{\alpha\hat\rho_k}{alpha * rho_k} if it is over-allocated (\eqn{N_k/m > \hat\rho_k}{N_k/m > rho_k}), \eqn{\hat\rho_k}{rho_k} if it is on target, and \eqn{\hat\rho_k + (1-\alpha)\sum_{j \in S}\hat\rho_j/|T|}{rho_k + (1 - alpha) sum_{j in S} rho_j / |T|} if it is under-allocated, where \eqn{S} and \eqn{T} are the sets of over- and under-allocated arms, \eqn{N_k/m} is the current allocation proportion and \eqn{\hat\rho}{rho} is the estimated target (Alkhnefr, Hu and Zhai, 2025, eq. (1)). For two arms this is the ERADE of Hu, Zhang and He (2009). ERADE attains the lower bound of the asymptotic variance of the allocation proportions, so its allocation is less variable than that of the DBCD.

Sequential monitoring. With \code{monitor = sqMonitor(t, spend)} (two arms only) the trial is analysed at the information times \code{t}. Looks are taken at the end of the first group in which at least \code{ceiling(t * ssn)} patients are enrolled. If one group passes several planned looks, only the last of them is analysed and the alpha of the skipped looks is not spent, which makes the test slightly conservative. At each look the Wald statistic \eqn{Z = (\hat\theta_1 - \hat\theta_2)/\sqrt{\hat v_1/n_1 + \hat v_2/n_2}}{Z = (theta_1 - theta_2)/sqrt(v_1/n_1 + v_2/n_2)} is computed from all patients enrolled so far, with the sample variances as variance estimates (missing responses are excluded), and the trial stops and rejects the null hypothesis as soon as \eqn{|Z|} reaches the boundary of \code{\link{sqBoundary}}. Zhu and Hu (2010) showed that under the DBCD the sequential statistics are asymptotically a Brownian motion in the information time, so alpha spending boundaries keep the type I error. Their theory covers two arms and the DBCD, and the same boundaries are used with ERADE. The result then also contains the stopping probability at each look and the expected sample size.
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
\item{data: allocation}{An \code{nsim} by \code{ssn} matrix with the treatment assignments of every simulated trial (positions after an early stop are \code{NA}).}
\item{type I error}{Only if \code{typeI = TRUE}. The rejection rate under the null hypothesis.}
\item{boundary}{Only if \code{monitor} is given. The boundaries for \eqn{|Z|} at the looks.}
\item{stopping probability}{Only if \code{monitor} is given. The proportion of simulated trials that stop at each look (the last look is the final analysis).}
\item{expected sample size}{Only if \code{monitor} is given. The mean number of enrolled patients.}
\item{data: stage}{Only if \code{monitor} is given. The look at which each simulated trial stopped.}
\item{data: sample size}{Only if \code{monitor} is given. The number of patients enrolled in each simulated trial.}
}
\references{
Alkhnefr, N., Hu, F. and Zhai, G. (2025). \emph{Efficient randomized adaptive designs for multi-arm clinical trials.} Statistical Methods in Medical Research, 34(9), 1886-1898. \doi{10.1177/09622802251362644}

Hu, F. and Zhang, L. X. (2004). \emph{Asymptotic properties of doubly adaptive biased coin designs for multitreatment clinical trials.} The Annals of Statistics, 32(1), 268-301. \doi{10.1214/aos/1079120137}

Hu, F., Zhang, L. X. and He, X. (2009). \emph{Efficient randomized-adaptive designs.} The Annals of Statistics, 37(5A), 2543-2560. \doi{10.1214/08-AOS655}

Tymofyeyev, Y., Rosenberger, W. F. and Hu, F. (2007). \emph{Implementing optimal allocation in sequential binary response experiments.} Journal of the American Statistical Association, 102(477), 224-234. \doi{10.1198/016214506000000906}

Zhai, G., Li, Y., Zhang, L. and Hu, F. (2024). \emph{Group response-adaptive randomization with delayed and missing responses.} Statistics in Medicine, 43(27), 5047-5059. \doi{10.1002/sim.10220}

Zhang, L. and Rosenberger, W. F. (2006). \emph{Response-adaptive randomization for clinical trials with continuous outcomes.} Biometrics, 62(2), 562-569. \doi{10.1111/j.1541-0420.2005.00496.x}

Zhu, H. and Hu, F. (2010). \emph{Sequential monitoring of response-adaptive randomized clinical trials.} The Annals of Statistics, 38(4), 2218-2241. \doi{10.1214/10-AOS796}
}

\examples{
theta = c(13, 4.0^2, 15, 2.5^2)
k = 2
gsize.param = 5
ssn = 120
res <- Group.DBCD_Cont(n0 = 20, theta = theta, k = k, gsize.param = gsize.param,
                       ssn = ssn, target.alloc = "Neyman", r = 2, nsim = 100,
                       mRate = NULL, alpha = 0.05)
res
}

\keyword{design}
\concept{group doubly adaptive biased coin design}
\concept{group response-adaptive randomization}
\concept{group doubly adaptive biased coin design with continuous response}
