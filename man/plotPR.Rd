\name{plotPR}
\alias{plotPR}
\title{
  Plot PR Curve
}
\description{
  A function that automatically produce a plot of the PR (Precision/Recall) Curve.
}
\usage{
plotPR(piNet, ...)
}
\arguments{
  \item{piNet}{
    the matrix generated from the \code{performanceIndex} function.
  }
  \item{\dots}{
    other parameters to be added to the \code{plot} function.
  }
}
\author{
  Luciano Garofano \email{lucianogarofano88@gmail.com}, Stefano Maria Pagnotta, Michele Ceccarelli
}
\seealso{
  \code{\link[synRNASeqNet:plotROC]{plotROC}}
}
\examples{
simData <- simulatedData(p = 10, n = 20, mu = 100, sigma = 0.25,
                        ppower = 0.73, noise = FALSE)
counts <- simData$counts
adjMat <- simData$adjMat

miML <- parMIEstimate(counts, method = "ML", unit = "nat", nchips = 2)
valML <- performanceIndex(miML, adjMat)

plotPR(valML, col = "red")
}
\keyword{plotPR}