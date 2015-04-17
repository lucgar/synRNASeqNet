\name{performanceNET}
\alias{performanceNET}
\title{
  Evalutate Performance Indices
}
\description{
  A function that computes the peformance indices (i.e. "Recall", "FPR", "Precision",
  "Accuracy", "Fscore") between a predicted adjacency matrix (\code{testNet})
  and a true one (\code{gsNet}) when the false interactions are already removed.
}
\usage{
performanceNET(predAdjMat, valAdjMat)
}
\arguments{
  \item{predAdjMat}{
    a predicted adjacency matrix to evaluate when the false interactions are already removed.
  }
  \item{valAdjMat}{
    the true adjacency matrix.
  }
}
\value{
  The \code{performanceIndex} function returns a data.frame with the considered
  performance indices.
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
  \link[synRNASeqNet:performanceIndex]{performanceIndex}}
}
\examples{
simData <- simulatedData(p = 10, n = 20, mu = 100, sigma = 0.25,
                        ppower = 0.73, noise = FALSE)
counts <- simData$counts
adjMat <- simData$adjMat

miML <- parMIEstimate(counts, method = "ML", unit = "nat", nchips = 2)
miML <- aracne.a(miML, eps = 0.5)

valML <- performanceNET(miML, adjMat)
}
\keyword{performanceNET}