\name{entropyML}
\alias{entropyML}
\title{
  Maximum Likelihood Entropy Estimate
}
\description{
  Computing the Maximul Likelihood Entropy Estimate of \code{cellCounts}.
}
\usage{
entropyML(cellCounts, unit = unit)
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
  The \code{entropyML} function returns the value of the entropy of that gene
  \code{H(X)} (or pair of genes \code{H(X,Y)}).
}
\author{
  Luciano Garofano \email{lucianogarofano88@gmail.com}, Stefano Maria Pagnotta, Michele Ceccarelli
}
\references{
  Paniski L. (2003). Estimation of Entropy and Mutual Information. \emph{Neural
  Computation}, vol. 15 no. 6 pp. 1191-1253.
  
  Meyer P.E., Laffitte F., Bontempi G. (2008). minet: A R/Bioconductor Package
  for Inferring Large Transcriptional Networks Using Mutual Information.
  \emph{BMC Bioinformatics} 9:461.
  
  Antos A., Kontoyiannis I. (2001). Convergence properties of functional estimates
  for discrete distributions. \emph{Random Structures and Algorithms}, vol. 19
  pp. 163-193.
  
  Strong S., Koberle R., de Ruyter van Steveninck R.R., Bialek W. (1998). Entropy
  and Information in Neural Spike Trains. \emph{Physical Review Letters}, vol.
  80 pp. 197-202.
}
\seealso{
  \code{\link[synRNASeqNet:entropyMM]{entropyMM}, \link[synRNASeqNet:entropyBayes]{entropyBayes},
  \link[synRNASeqNet:entropyCS]{entropyCS}, \link[synRNASeqNet:entropyShrink]{entropyShrink}}
}
\examples{
simData <- simulatedData(p = 50, n = 100, mu = 100, sigma = 0.25,
                        ppower = 0.73, noise = FALSE)
cellCounts <- table(simData$counts[1, ])
eML <- entropyML(cellCounts, unit = "nat")
}
\keyword{entropyML}
\keyword{thetaML}