\name{performanceIndex}
\alias{performanceIndex}
\title{
  Evalutate Performance Indices
}
\description{
  A function that computes the peformance indices (i.e. "Recall", "FPR", "Precision",
  "Accuracy", "Fscore") between a predicted adjacency matrix (\code{testNet})
  and a true one (\code{gsNet}) for each mutual information threshold.
}
\usage{
performanceIndex(testNet, gsNet)
}
\arguments{
  \item{testNet}{
    a predicted adjacency matrix to evaluate.
  }
  \item{gsNet}{
    the true adjacency matrix.
  }
}
\value{
  The \code{performanceIndex} function returns a data.frame with the considered
  performance indices for each mutual information threshold.
}
\author{
  Luciano Garofano \email{lucianogarofano88@gmail.com}, Stefano Maria Pagnotta, Michele Ceccarelli
}
\references{
  Stehman, S.V. (1997). Selecting and interpreting measures of thematic classification
  accuracy. \emph{Remote Sensing of Environment} 62 (1): 77-89.
}
\seealso{
  \code{\link[synRNASeqNet:simulatedData]{simulatedData},
  \link[synRNASeqNet:mainNetFunction]{mainNetFunction},
  \link[synRNASeqNet:performanceNET]{performanceNET}}
}
\examples{
simData <- simulatedData(p = 10, n = 20, mu = 100, sigma = 0.25,
                        ppower = 0.73, noise = FALSE)
counts <- simData$counts
adjMat <- simData$adjMat

miML <- parMIEstimate(counts, method = "ML", unit = "nat", nchips = 2)
valML <- performanceIndex(miML, adjMat)
}
\keyword{performanceIndex}