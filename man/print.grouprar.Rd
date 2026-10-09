\name{print.grouprar}
\alias{print.grouprar}
\alias{summary.grouprar}
\alias{print.summary.grouprar}
\title{
Print and Summarize Simulation Results
}
\description{
Methods for the objects of class \code{"grouprar"} returned by the simulation functions of the package, such as \code{\link{DBCD_Bin}} or \code{\link{WeiUrn}}.
}
\usage{
\method{print}{grouprar}(x, \dots)

\method{summary}{grouprar}(object, \dots)

\method{print}{summary.grouprar}(x, digits = 3, \dots)
}
\arguments{
  \item{x}{An object of class \code{"grouprar"} (for \code{print}) or \code{"summary.grouprar"} (for \code{print.summary.grouprar}).}
  \item{object}{An object of class \code{"grouprar"}.}
  \item{digits}{The number of decimal places to print. The default is 3.}
  \item{\dots}{Further arguments passed to or from other methods.}
}
\details{
\code{print} shows the summary of the object. The full simulation results, including the per-simulation data, remain available as list elements, for example \code{x[["data: allocation"]]}.
}
\value{
\code{summary} returns an object of class \code{"summary.grouprar"}, a list with the elements
\item{method}{The name of the procedure.}
\item{sample size}{The sample size.}
\item{nsim}{The number of simulated trials.}
\item{arms}{A data frame with one row per arm: the true parameters, the mean allocation proportion and its standard deviation over the simulations.}
\item{stats}{A named vector with the failure rate (or mean response for continuous designs) and its standard deviation, the power, and, when available, the type I error, the trial and enrollment durations, and the expected sample size.}
\item{stopping probability}{For monitored designs, the proportion of simulated trials that stop at each look.}

The \code{print} methods return their argument invisibly.
}
\seealso{
\code{\link{DBCD_Bin}}, \code{\link{CRDesign}}.
}
\examples{
res <- CRDesign(k = 3, p = c(0.6, 0.7, 0.8), ssn = 100, nsim = 50, seed = 1)
res
s <- summary(res)
s$arms
print(s, digits = 2)
}
\keyword{print}
