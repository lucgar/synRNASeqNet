\name{synRNASeqNet-package}
\alias{synRNASeqNet-package}
\alias{synRNASeqNet}
\docType{package}
\title{
  Synthetic RNA-Seq Network Generation and Mutual Information Estimates
}
\description{
  It implements various estimators of mutual information, such as
  the maximum likelihood and the Millow-Madow estimator, various Bayesian
	estimators, the shrinkage estimator, and the Chao-Shen estimator. It also
  offers wrappers to the kNN and kernel density estimators.
  
  Furthermore, it provides various index of performance evaluation such as
  precision, recall, FPR, F-Score, ROC-PR Curves and so on.
  
  Lastly, it provides a brand new way of generating synthetic RNA-Seq Network
  with known dependence structure.
}

\section{Index}{
  
  \tabular{ll}{
    \link[synRNASeqNet:synRNASeqNet-package]{synRNASeqNet-package}: \tab Synthetic RNA-Seq Network Generation and Mutual Information Estimates\cr
    \link[synRNASeqNet:parMIEstimate]{parMIEstimate}: \tab Parallel Mutual Information Estimation\cr
    \link[synRNASeqNet:parEntropyEstimate]{parEntropyEstimate}: \tab Parallel Entropy Estimation\cr
    \link[synRNASeqNet:entropyML]{entropyML}: \tab Maximum Likelihood Entropy Estimate\cr
    \link[synRNASeqNet:entropyMM]{entropyMM}: \tab Miller-Madow corrected Entropy Estimate\cr
    \link[synRNASeqNet:entropyBayes]{entropyBayes}: \tab Bayesian Entropy Estimate\cr
    \link[synRNASeqNet:entropyCS]{entropyCS}: \tab Chao-Shen Entropy Estimate\cr
    \link[synRNASeqNet:entropyShrink]{entropyShrink}: \tab James-Stein Shrinkage Entropy Estimate\cr
    \link[synRNASeqNet:parMIKD]{parMIKD}: \tab Parallel Kernel Density Mutual Information Estimate\cr
    \link[synRNASeqNet:simulatedData]{simulatedData}: \tab Random Generation Networks for RNA-Seq Data\cr
    \link[synRNASeqNet:mainNetFunction]{mainNetFunction}: \tab Main Estimation and Evaluation Function\cr
    \link[synRNASeqNet:plotROC]{plotROC}: \tab Plot ROC Curve\cr
    \link[synRNASeqNet:plotPR]{plotPR}: \tab Plot PR Curve\cr
    \link[synRNASeqNet:aucDisc]{aucDisc}: \tab Calculate Area Under a (ROC/PR) Curve\cr
    \link[synRNASeqNet:performanceIndex]{performanceIndex}: \tab Evalutate Performance Indices\cr
    \link[synRNASeqNet:performanceNET]{performanceNET}: \tab Evalutate Performance Indices\cr
    \link[synRNASeqNet:YoudenIndex]{YoudenIndex}: \tab Youden's Index\cr
    \link[synRNASeqNet:Likelihoods]{Likelihoods}: \tab Likelihood Indices\cr
    \link[synRNASeqNet:DiscriminantPower]{DiscriminantPower}: \tab Discriminant Power\cr
  }
}

\author{
  Luciano Garofano \email{lucianogarofano88@gmail.com}, Stefano Maria Pagnotta, Michele Ceccarelli
}
\keyword{package}