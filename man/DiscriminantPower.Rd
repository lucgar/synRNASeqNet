\name{DiscriminantPower}
\alias{DiscriminantPower}
\title{
  Discriminant Power
}
\description{
  The Discriminant Power summarizes sensitivity and specificity.
}
\usage{
DiscriminantPower(resTable)
}
\arguments{
  \item{resTable}{
    the data.frame generated from the \code{performanceNET} (or one of the row
    of the matrix generated from the \code{performanceIndex}) function.
  }
}
\value{
  The \code{DiscriminantPower} function returns the value of the Discriminant Power
  for a given predicted adjacency matrix.
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
  \code{\link[synRNASeqNet:YoudenIndex]{YoudenIndex}, \link[synRNASeqNet:Likelihoods]{Likelihoods}}
}
\examples{
simData <- simulatedData(p = 10, n = 20, mu = 100, sigma = 0.25,
                        ppower = 0.73, noise = FALSE)
counts <- simData$counts
adjMat <- simData$adjMat

miML <- parMIEstimate(counts, method = "ML", unit = "nat", nchips = 2)
miML <- aracne.a(miML, eps = 0.5)

valML <- performanceNET(miML, adjMat)
DP <- DiscriminantPower(valML)
}
\keyword{DiscriminantPower}