\name{parEntropyEstimate}
\alias{parEntropyEstimate}
\title{
  Parallel Entropy Estimation
}
\description{
  A function that computes the entropy between all pairs of rows (or
  specified ones) of matrix \code{counts} using the indirect methods.
}
\usage{
parEntropyEstimate(idx, method = method, unit = unit,
                  priorHyperParam = priorHyperParam,
                   shrinkageTarget = shrinkageTarget, boot = boot)
}
\arguments{
  \item{idx}{
    the index of the cell which corresponds to the interaction going to be esimated.
  }
  \item{method}{
    a character string indicating which estimate is to be computed.
    One of \code{"ML"} (Maximum Likelihood Estimator, default), \code{"MM"} (Miller-Madow
    corrected Estimator), \code{"Bayes"} (Bayesian Estimators), \code{"CS"}
    (Chao-Shen Estimator), \code{"Shrink"} (James-Stein shrinkage Estimator),
    \code{"KD"} (kernel Density Estimator), or \code{"KNN"} (k-Nearest Neighbor
    Estimator), can be abbreviated. For the \code{"Bayes"} estimate it is needed
    to specify also which \code{priorHyperParam} is to be used; for \code{"Shrink"}
    is optional to specify values for the \code{shrinkageTarget} parameter; for
    \code{"KNN"} is needed to specify also the number of nearest neighbors \code{k}.
  }
  \item{unit}{
    the unit in which mutual information is measured. One of \code{"bit"}
    (log2, default), \code{"ban"} (log10) or \code{"nat"} (natural units).
  }
  \item{priorHyperParam}{
    the prior distribution type for the Bayes estimation. One of \code{"Jeffreys"}
    (default, Jeffreys Prior, Krichevsky and Trofimov Estimator), \code{"BLUnif"}
    (Bayes-Laplace uniform Prior, Holste Estimator), \code{"Perks"} (Perks Prior,
    Schurmann and Grassberger Estimator), or \code{"MiniMax"} (MiniMax Prior), can be
    abbreviated.
  }
  \item{shrinkageTarget}{
    shrinkage target frequencies. If not specified (default) it is estimated in
    a James-Stein-type fashion (uniform distribution).
  }
  \item{boot}{
    logical (\code{FALSE} as default). Used for calculating a null distribution
    in order to evaluate if such a interaction is true or obtained by chance.
  }
}
\details{
  Internal of \code{parMIEstimate}.
}
\value{
  The \code{parEntropyEstimate} function returns the value of the entropy of
  that pair of genes \code{H(X,Y)}.
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
  
  Miller G.A. (1955). Note on the bias of information estimates. \emph{Information
  Theory in Psychology}, II-B pp. 95-100.
  
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
  
  Chao A., Shen T.J. (2003). Nonparametric estimation of Shannon's index diversity
  when there are unseen species. \emph{Environmental and Ecological Statistics},
  vol. 10 pp. 429-443.
  
  James W., Stein C. (1961). Estimation with Quadratic Loss. \emph{Proceedings
  of the Fourth Berkeley Symposium on Mathematical Statistics and Probability},
  vol. 1 pp. 361-379.
}
\seealso{
  \code{\link[synRNASeqNet:parMIEstimate]{parMIEstimate}, \link[synRNASeqNet:parMIKD]{parMIKD}}
}
\keyword{parEntropyEstimate}