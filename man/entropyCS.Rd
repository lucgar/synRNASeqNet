\name{entropyCS}
\alias{entropyCS}
\title{
  Chao-Shen Entropy Estimate
}
\description{
  Computing the Chao-Shen Entropy Estimate of \code{cellCounts}.
}
\usage{
entropyCS(cellCounts, unit = unit)
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
}
\value{
  The \code{entropyCS} function returns the value of the entropy of that gene
  \code{H(X)} (or pair of genes \code{H(X,Y)}).
}
\author{
  Luciano Garofano \email{lucianogarofano88@gmail.com}, Stefano Maria Pagnotta, Michele Ceccarelli
}
\references{
  Chao A., Shen T.J. (2003). Nonparametric estimation of Shannon's index diversity
  when there are unseen species. \emph{Environmental and Ecological Statistics},
  vol. 10 pp. 429-443.
}
\seealso{
  \code{\link[synRNASeqNet:entropyML]{entropyML}, \link[synRNASeqNet:entropyMM]{entropyMM},
  \link[synRNASeqNet:entropyBayes]{entropyBayes}, \link[synRNASeqNet:entropyShrink]{entropyShrink}}
}
\examples{
simData <- simulatedData(p = 50, n = 100, mu = 100, sigma = 0.25,
                        ppower = 0.73, noise = FALSE)
cellCounts <- table(simData$counts[1, ])
eCS <- entropyCS(cellCounts, unit = "nat")
}
\keyword{entropyCS}
\keyword{thetaGT}