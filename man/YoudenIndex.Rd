\name{YoudenIndex}
\alias{YoudenIndex}
\title{
  Youden's Index
}
\description{
  The Youden's Index evaluates the algorithm's ability to avoid failure.
}
\usage{
YoudenIndex(resTable)
}
\arguments{
  \item{resTable}{
    the data.frame generated from the \code{performanceNET} (or one of the row
    of the matrix generated from the \code{performanceIndex}) function.
  }
}
\value{
  The \code{YoudenIndex} function returns the value of the Youden's Index for a
  given predicted adjacency matrix.
}
\author{
  Luciano Garofano \email{lucianogarofano88@gmail.com}, Stefano Maria Pagnotta, Michele Ceccarelli
}
\references{
  Sokolova M., Japkowicz N., Szpakowicz S. (2006). Beyond Accuracy, F-Score and
  ROC: A Family of Discriminant Measures for Performance Evaluation. \emph{Lecture
  Notes in Computer Science} Volume 4304, 2006, pp 1015-1021.
}
\seealso{
  \code{\link[synRNASeqNet:Likelihoods]{Likelihoods},
  \link[synRNASeqNet:DiscriminantPower]{DiscriminantPower}}
}
\examples{
simData <- simulatedData(p = 10, n = 20, mu = 100, sigma = 0.25,
                        ppower = 0.73, noise = FALSE)
counts <- simData$counts
adjMat <- simData$adjMat

miML <- parMIEstimate(counts, method = "ML", unit = "nat", nchips = 2)
miML <- aracne.a(miML, eps = 0.5)

valML <- performanceNET(miML, adjMat)
YI <- YoudenIndex(valML)
}
\keyword{YoudenIndex}