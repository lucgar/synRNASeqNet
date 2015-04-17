\name{simulatedData}
\alias{simulatedData}
\title{
  Random Generation Networks for RNA-Seq Data
}
\description{
  A function that use a stochastic BA-modelfor building a graph and
  the simulated RNA-Seq counts (from a Poisson (multivariate or over-dispersed)
  distribution) that encode the underlying graph structure.
}
\usage{
simulatedData(p = 50, n = 100, mu = 100, sigma = 0.25, ppower = 1,
              noise = F, seed = NULL)
}
\arguments{
  \item{p}{
    the number of genes in the networks.
  }
  \item{n}{
    the number of simulated samples for each gene.
  }
  \item{mu}{
    the average mean of the simulated Poisson distributions.
  }
  \item{sigma}{
    the over-dispersed sd value in the case of over-dispersed Poisson simulation.
  }
  \item{ppower}{
    the power of the preferential attachment for the BA-model.
  }
  \item{noise}{
    logical. Should same noise be added to the data or not?
  }
  \item{seed}{
    a single value, interpreted as an integer, in order to control the simulated data.
  }
}
\value{
  \item{graph}{
    the graph generated with the BA-model.
  }
  \item{adjMat}{
    the related adjacency matrix that encodes the underlying graph structure.
  }
  \item{counts}{
    the simulated RNA-Seq counts matrix.
  }
}
\author{
  Luciano Garofano \email{lucianogarofano88@gmail.com}, Stefano Maria Pagnotta, Michele Ceccarelli
}
\references{
  Barabasi A.L., Albert R. (1999). Emergence of scaling in random networks.
  \emph{Science}, 286 509-512.
  
  Gallopin M., Rau A., Jaffrezic F. (2013). A Hierarchical Poisson Log-Normal
  Model for Network Inference from RNA Sequencing Data. \emph{PLOSone}.
}
\seealso{
  \code{\link[synRNASeqNet:mainNetFunction]{mainNetFunction}}
}
\examples{
simData <- simulatedData(p = 50, n = 100, mu = 100, sigma = 0.25,
                        ppower = 0.73, noise = FALSE)
plot(simData$graph)

counts <- simData$counts
adjMat <- simData$adjMat
}
\keyword{simulatedData}