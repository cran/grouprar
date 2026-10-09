\name{grouprar}
\alias{grouprar}
\alias{grouprar-package}
\docType{package}
\title{
  grouprar-package: Group Response-Adaptive Randomization for Clinical Trials
}
\description{
  Implements group response-adaptive randomization procedures, which include standard (non-group) response-adaptive randomization methods as special cases. The package also handles delayed and missing responses, which broadens its use in real-world trials. It offers functions for simulating a variety of response-adaptive randomization procedures, to help guide the choice of design for a clinical trial.

  Urn designs: \code{\link{RPWRule}}, \code{\link{WeiUrn}}, \code{\link{PolyaUrn}}, \code{\link{DLRule}}, \code{\link{GDLRule}}, \code{\link{Bai.Hu.Shen.Urn}} and \code{\link{BirthDeathUrn}}, with \code{\link{CRDesign}} (complete randomization) as a non-adaptive benchmark.

  Doubly adaptive biased coin designs (DBCD) for binary and continuous responses, with immediate or delayed responses, for individual patients or groups of patients: \code{\link{DBCD_Bin}}, \code{\link{DBCD_Cont}}, \code{\link{dyldDBCD_Bin}}, \code{\link{dyldDBCD_Cont}}, \code{\link{Group.DBCD_Bin}}, \code{\link{Group.DBCD_Cont}}, \code{\link{Group.dyldDBCD_Bin}} and \code{\link{Group.dyldDBCD_Cont}}. These designs support any number of arms, the allocation function of Hu and Zhang (2004) or the multi-arm efficient randomized adaptive design (ERADE), several target allocations including the k-arm optimal allocations of Tymofyeyev, Rosenberger and Hu (2007) with a lower bound, and missing responses.

  Further tools: \code{\link{nextAlloc}} computes the allocation probabilities of the next patient or group in an ongoing trial, and \code{\link{sqMonitor}} and \code{\link{sqBoundary}} add group sequential monitoring with alpha spending to two-arm designs (Zhu and Hu, 2010). Every design accepts a user-supplied test (\code{test.fun}), can estimate the type I error (\code{typeI}) and takes a \code{seed}. The results are objects of class \code{"grouprar"} with \code{print} and \code{summary} methods (\code{\link{print.grouprar}}), and contain the allocation sequences of all simulated trials and, for delayed designs, the trial duration.

  For details of the methods and algorithms, see the following references:
    Wei, L. J. (1979) \doi{10.1214/aos/1176344614};
    Wei, L. J. and Durham, S. (1978) \doi{10.1080/01621459.1978.10480109};
    Durham, S. D., Flournoy, N. and Li, W. (1998) \doi{10.2307/3315771};
    Ivanova, A., Rosenberger, W. F., Durham, S. D. and Flournoy, N. (2000) \url{https://www.jstor.org/stable/25053121};
    Bai, Z. D., Hu, F. and Shen, L. (2002) \doi{10.1006/jmva.2001.1987};
    Ivanova, A. (2003) \doi{10.1007/s001840200220};
    Hu, F. and Zhang, L. X. (2004) \doi{10.1214/aos/1079120137};
    Hu, F. and Rosenberger, W. F. (2006, ISBN:978-0-471-65396-7);
    Zhang, L. X., Chan, W. S., Cheung, S. H. and Hu, F. (2007) \url{https://www.jstor.org/stable/26432528};
    Zhang, L. and Rosenberger, W. F. (2006) \doi{10.1111/j.1541-0420.2005.00496.x};
    Tymofyeyev, Y., Rosenberger, W. F. and Hu, F. (2007) \doi{10.1198/016214506000000906};
    Hu, F., Zhang, L. X., Cheung, S. H. and Chan, W. S. (2008) \doi{10.1002/cjs.5550360404};
    Hu, F., Zhang, L. X. and He, X. (2009) \doi{10.1214/08-AOS655};
    Zhu, H. and Hu, F. (2010) \doi{10.1214/10-AOS796};
    Zhai, G., Li, Y., Zhang, L. and Hu, F. (2024) \doi{10.1002/sim.10220};
    Alkhnefr, N., Hu, F. and Zhai, G. (2025) \doi{10.1177/09622802251362644}.
}

\references{
Alkhnefr, N., Hu, F. and Zhai, G. (2025). \emph{Efficient randomized adaptive designs for multi-arm clinical trials.} Statistical Methods in Medical Research, 34(9), 1886-1898.

Bai, Z. D., Hu, F. and Shen, L. (2002). \emph{An adaptive design for multi-arm clinical trials.} Journal of Multivariate Analysis, 81(1), 1-18.

Durham, S. D., Flournoy, N. and Li, W. (1998). \emph{A sequential design for maximizing the probability of a favourable response.} Canadian Journal of Statistics, 26(3), 479-495.

Hu, F. and Rosenberger, W. F. (2003). \emph{Optimality, variability, power: evaluating response-adaptive randomization procedures for treatment comparisons.} Journal of the American Statistical Association, 98(463), 671-678.

Hu, F. and Rosenberger, W. F. (2006). \emph{The Theory of Response-Adaptive Randomization in Clinical Trials.} John Wiley & Sons.

Hu, F. and Zhang, L. X. (2004). \emph{Asymptotic properties of doubly adaptive biased coin designs for multitreatment clinical trials.} The Annals of Statistics, 32(1), 268-301.

Hu, F., Zhang, L. X. and He, X. (2009). \emph{Efficient randomized-adaptive designs.} The Annals of Statistics, 37(5A), 2543-2560.

Hu, F., Zhang, L. X., Cheung, S. H. and Chan, W. S. (2008). \emph{Doubly adaptive biased coin designs with delayed responses.} Canadian Journal of Statistics, 36(4), 541-559.

Ivanova, A. (2003). \emph{A play-the-winner-type urn design with reduced variability.} Metrika, 58, 1-13.

Ivanova, A., Rosenberger, W. F., Durham, S. D. and Flournoy, N. (2000). \emph{A birth and death urn for randomized clinical trials: asymptotic methods.} Sankhya, Series B, 62(1), 104-118.

Tymofyeyev, Y., Rosenberger, W. F. and Hu, F. (2007). \emph{Implementing optimal allocation in sequential binary response experiments.} Journal of the American Statistical Association, 102(477), 224-234.

Wei, L. J. (1979). \emph{The generalized Polya's urn design for sequential medical trials.} The Annals of Statistics, 7(2), 291-296.

Wei, L. J. and Durham, S. (1978). \emph{The randomized play-the-winner rule in medical trials.} Journal of the American Statistical Association, 73(364), 840-843.

Zelen, M. (1969). \emph{Play the winner rule and the controlled clinical trial.} Journal of the American Statistical Association, 64(325), 131-146.

Zhai, G., Li, Y., Zhang, L. and Hu, F. (2024). \emph{Group response-adaptive randomization with delayed and missing responses.} Statistics in Medicine, 43(27), 5047-5059.

Zhang, L. and Rosenberger, W. F. (2006). \emph{Response-adaptive randomization for clinical trials with continuous outcomes.} Biometrics, 62(2), 562-569.

Zhang, L. X., Chan, W. S., Cheung, S. H. and Hu, F. (2007). \emph{A generalized drop-the-loser urn for clinical trials with delayed responses.} Statistica Sinica, 17(1), 387-409.

Zhu, H. and Hu, F. (2010). \emph{Sequential monitoring of response-adaptive randomized clinical trials.} The Annals of Statistics, 38(4), 2218-2241.
}

\author{
Guannan Zhai \email{guannanzhai1996@gmail.com}; Feifang Hu \email{feifang@gwu.edu}.
}
\keyword{package}
