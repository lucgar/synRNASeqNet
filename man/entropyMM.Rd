\name{entropyMM}
\alias{entropyMM}
\title{
  Miller-Madow corrected Entropy Estimate
}
\description{
  Computing the Miller-Madow corrected Entropy Estimate of \code{cellCounts}.
}
\usage{
entropyMM(cellCounts, unit = unit)
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
  The \code{entropyMM} function returns the value of the entropy of that gene
  \code{H(X)} (or pair of genes \code{H(X,Y)}).
}
\author{
  Luciano Garofano \email{lucianogarofano88@gmail.com}, Stefano Maria Pagnotta, Michele Ceccarelli
}
\references{
  Miller G.A. (1955). Note on the bias of information estimates. \emph{Information
  Theory in Psychology}, II-B pp. 95-100.
}
\seealso{
  \code{\link[synRNASeqNet:entropyMM]{entropyML}, \link[synRNASeqNet:entropyMM]{entropyBayes},
  \link[synRNASeqNet:entropyMM]{entropyCS}, \link[synRNASeqNet:entropyMM]{entropyShrink}}
}
\examples{
simData <- simulatedData(p = 50, n = 100, mu = 100, sigma = 0.25,
                        ppower = 0.73, noise = FALSE)
cellCounts <- table(simData$counts[1, ])
eMM <- entropyMM(cellCounts, unit = "nat")
}
\keyword{entropyMM}