\name{entropyShrink}
\alias{entropyShrink}
\title{
  James-Stein Shrinkage Entropy Estimate
}
\description{
  Computing the James-Stein Shrinkage Entropy Estimate of \code{cellCounts}.
}
\usage{
entropyShrink(cellCounts, unit = unit, shrinkageTarget = shrinkageTarget)
}
\arguments{
  \item{cellCounts}{
    an integer vector (or matrix) representing the number of times each particular
    count is obtained.
  }
  \item{unit}{
    the unit in which entropy is measured. One of \code{"bit"} (log2, default),
    \code{"ban"} (log10) or \code{"nat"} (natural units).
  }
  \item{shrinkageTarget}{
    shrinkage target frequencies. If not specified (default) it is estimated in
    a James-Stein-type fashion (uniform distribution).
  }
}
\value{
  The \code{entropyShrink} function returns the value of the entropy of that gene
  \code{H(X)} (or pair of genes \code{H(X,Y)}).
}
\author{
  Luciano Garofano \email{lucianogarofano88@gmail.com}, Stefano Maria Pagnotta, Michele Ceccarelli
}
\references{
  James W., Stein C. (1961). Estimation with Quadratic Loss. \emph{Proceedings
  of the Fourth Berkeley Symposium on Mathematical Statistics and Probability},
  vol. 1 pp. 361-379.
}
\seealso{
  \code{\link[synRNASeqNet:entropyShrink]{entropyML},
  \link[synRNASeqNet:entropyShrink]{entropyMM}, \link[synRNASeqNet:entropyShrink]{entropyBayes},
  \link[synRNASeqNet:entropyShrink]{entropyCS}}
}
\examples{
simData <- simulatedData(p = 50, n = 100, mu = 100, sigma = 0.25,
                        ppower = 0.73, noise = FALSE)
cellCounts <- table(simData$counts[1, ])
eShrink <- entropyShrink(cellCounts, unit = "nat", shrinkageTarget = NULL)
}
\keyword{entropyShrink}
\keyword{thetaShrink}
\keyword{shrinkageIntensity}