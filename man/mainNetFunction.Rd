\name{mainNetFunction}
\alias{mainNetFunction}
\title{
  Main Estimation and Evaluation Function
}
\description{
  A function that computes the mutual information between all pairs of rows (or
  specified ones) of matrix \code{counts} using all the 10 different estimation
  methods and evalute their performances.
}
\usage{
mainNetFunction(counts, adjMat, nchips, plotPath = "",
                tfList = NULL)
}
\arguments{
  \item{counts}{
    a numeric matrix (for the reconstruction of gene regulatory networks, genes
    on rows and samples on columns).
  }
  \item{adjMat}{
    the adjacency matrix that encodes the graph structure that is going to be
    predicted.
  }
  \item{nchips}{
    the number of cpu's to be used for making the parallel calculation.
  }
  \item{plotPath}{
    the folder in which the plot will be saved.
  }
  \item{tfList}{
    the character vector specifying which genes from the rownames of the \code{counts}
    matrix is to be used as transcription factor for network reconstruction.
  }
}
\value{
  \item{miEst}{
    a list containing the estimates of all methods.
  }
  \item{valMet}{
    a list contatining the performance indices (i.e. "Recall", "FPR", "Precision",
    "Accuracy", "Fscore") calculated in all methods and usable for creating
    curves like ROC and PR.
  }
  \item{resTable}{
    a matrix with the best performces for each method.
  }
}
\author{
  Luciano Garofano \email{lucianogarofano88@gmail.com}, Stefano Maria Pagnotta, Michele Ceccarelli
}
\references{
  Stehman, S.V. (1997). Selecting and interpreting measures of thematic classification
  accuracy. \emph{Remote Sensing of Environment} 62 (1): 77-89.
}
\seealso{
  \code{\link[synRNASeqNet:simulatedData]{simulatedData}}
}
\examples{
simData <- simulatedData(p = 5, n = 10, mu = 50, sigma = 0.25,
                        ppower = 0.73, noise = FALSE)
counts <- simData$counts
adjMat <- simData$adjMat

#netData <- mainNetFunction(counts, adjMat, nchips = 2)
}
\keyword{mainNetFunction}