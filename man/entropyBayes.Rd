\name{entropyBayes}
\alias{entropyBayes}
\title{
  Bayesian Entropy Estimate
}
\description{
  Computing the Maximul Likelihood Entropy Estimate of \code{cellCounts}.
}
\usage{
entropyBayes(cellCounts, unit = unit,
            priorHyperParam = priorHyperParam)
}\arguments{
  \item{cellCounts}{
    an integer vector (or matrix) representing the number of times each particular
    count is obtained.
  }
  \item{unit}{
    the unit in which entropy is measured. One of \code{"bit"} (log2, default),
    \code{"ban"} (log10) or \code{"nat"} (natural units).
  }
  \item{priorHyperParam}{
    the prior distribution type for the Bayes estimation. One of \code{"Jeffreys"}
    (default, Jeffreys Prior, Krichevsky and Trofimov Estimator), \code{"BLUnif"}
    (Bayes-Laplace uniform Prior, Holste Estimator), \code{"Perks"} (Perks Prior,
    Schurmann and Grassberger Estimator), or \code{"MiniMax"} (MiniMax Prior), can be
    abbreviated.
  }
}
\value{
  The \code{entropyBayes} function returns the value of the entropy of that gene
  \code{H(X)} (or pair of genes \code{H(X,Y)}).
}
\author{
  Luciano Garofano \email{lucianogarofano88@gmail.com}, Stefano Maria Pagnotta, Michele Ceccarelli
}
\references{
  Jeffreys H. (1946). An invariant form for the prior probability in estimation
  problems. \emph{Proceedings of the Royal Society of London}, vol. 186 no. 1007
  pp. 453-461.
  
  Krichevsky R.E., Trofimov V.K. (1981). The performance of universal encoding.
  \emph{IEEE Transactions on Information Theory}, vol. 27 pp. 199-207.
  
  Holste D., Hertzel H. (1998). Bayes' estimators of generalized entropies.
  \emph{Journal of Physics A}, vol. 31 pp. 2551-2566.
  
  Perks W. (1947). Some observations on inverse probability including a new
  indifference rule. \emph{Journal of the Institute of Actuaries}, vol. 73 pp.
  285-334.
  
  Schurmann T., Grassberg P. (1996). Entropy estimation of symbol sequences.
  \emph{Chaos}, vol. 6 pp. 414-427.
  
  Trybula S. (1958). Some problems of simultaneous minimax estimation. \emph{The
  Annals of Mathematical Statistics}, vol. 29 pp. 245-253.
}
\seealso{
  \code{\link[synRNASeqNet:entropyML]{entropyML}, \link[synRNASeqNet:entropyMM]{entropyMM},
  \link[synRNASeqNet:entropyCS]{entropyCS}, \link[synRNASeqNet:entropyShrink]{entropyShrink}}
}
\examples{
simData <- simulatedData(p = 50, n = 100, mu = 100, sigma = 0.25,
                        ppower = 0.73, noise = FALSE)
cellCounts <- table(simData$counts[1, ])
eBJ <- entropyBayes(cellCounts, unit = "nat", priorHyperParam = "Jeffreys")
eBB <- entropyBayes(cellCounts, unit = "nat", priorHyperParam = "BLUnif")
eBP <- entropyBayes(cellCounts, unit = "nat", priorHyperParam = "Perks")
eBM <- entropyBayes(cellCounts, unit = "nat", priorHyperParam = "MiniMax")
}
\keyword{entropyBayes}
\keyword{thetaBayes}
\keyword{beta_k}