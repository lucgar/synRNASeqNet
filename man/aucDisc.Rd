\name{aucDisc}
\alias{aucDisc}
\title{
  Calculate Area Under a (ROC/PR) Curve
}
\description{
  aucDisc calculates the Area under a Curve giving the coordinate of each point.
}
\usage{
  aucDisc(fpr, tpr)
}
\arguments{
  \item{fpr}{
    the false positive rate (precision/abscissa) that forms the curve.
  }
  \item{tpr}{
    the true positive rate (recall/ordinate) that forms the curve.
  }
}
\details{
  It uses the 2nd order Taylor series in order to calculate the area under the
  curve.
}
\value{
  The \code{aucDisc} returns a real number representing the Area Under the given Curve.
}
\author{
  Luciano Garofano \email{lucianogarofano88@gmail.com}, Stefano Maria Pagnotta, Michele Ceccarelli
}
\references{
   Swets, J.A. (1996). Signal detection theory and ROC analysis in psychology
   and diagnostics: collected papers. \emph{Lawrence Erlbaum Associates},
   Mahwah, NJ.
}
\seealso{
  \code{\link[synRNASeqNet:performanceIndex]{performanceIndex}}
}
\examples{
simData <- simulatedData(p = 50, n = 100, mu = 100, sigma = 0.25,
                        ppower = 0.73, noise = FALSE)
counts <- simData$counts
adjMat <- simData$adjMat

miML <- parMIEstimate(counts, method = "ML", unit = "nat", nchips = 2)
valML <- performanceIndex(miML, adjMat)

AUROC <- aucDisc(valML[, "FPR"], valML[, "Recall"])
AUPR <- aucDisc(valML[, "Recall"], valML[, "Precision"])
}
\keyword{aucDisc}
%\keyword{ ~kwd2 }
