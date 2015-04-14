\name{parMIKD}
\alias{parMIKD}
\title{
  Parallel Kernel Density Mutual Information Estimate
}
\description{
  A function that computes the mutual information between all pairs of rows (or
  specified ones) of matrix \code{counts} using Kernel Density Estimate. It is
  a wrapper of the functions of the package \code{GenKern}.
}
\usage{
parMIKD(idx)
}
\arguments{
  \item{idx}{
    a numeric matrix (for the reconstruction of gene regulatory networks, genes
    on rows and samples on columns).
  }
}
\details{
  Internal of \code{parMIEstimate}.
}
\value{
  The \code{parMIKD} function returns the value of the mutual information of
  that pair of genes \code{I(X,Y)}.
}
\author{
  Luciano Garofano \email{lucianogarofano88@gmail.com}, Stefano Maria Pagnotta, Michele Ceccarelli
}
\references{
  Moon Y., Rajagopalan B., Lall U. (1995). Estimation of mutual information using
  kernel density estimators. \emph{Physical Review E}, vol. 52 n. 3 pp. 2318-2321.
}
\seealso{
  \code{\link[synRNASeqNet:parEntropyEstimate]{parEntropyEstimate},
  \link[synRNASeqNet:parMIEstimate]{parMIEstimate}}
}
\examples{
simData <- simulatedData(p = 10, n = 50, mu = 100, sigma = 0.25,
                        ppower = 0.73, noise = FALSE)
counts <- simData$counts
adjMat <- simData$adjMat

miKD <- parMIEstimate(counts, method = "KD", nchips = 2)
}
\keyword{parMIKD}