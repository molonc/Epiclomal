\name{REpiclomal-package}
\alias{REpiclomal-package}
\alias{REpiclomal}
\docType{package}
\title{
  R Package for Epiclomal
}
\description{
  R Package for Epiclomal. Includes non-probabilistic clustering methods to cluster methylation data, evaluation software for epiclomal clustering algorithm. Used by the epiclomal workflow.
}
\details{
  Clustering methods include EuclideanClust, HammingClust (PBAL), PearsonClust, DensityCut
}
\author{
Camila P.E. de Souza, Mirela Andronescu

Maintainer: Patricia Ye <pye@bccrc.ca>
}
\references{
  de Souza CPE, Andronescu M, Masud T, Kabeer F, Biele J, Laks E, et al. Epiclomal: probabilistic clustering of sparse single-cell DNA methylation data. bioRxiv. 2018;414482.
}
\seealso{
\code{\link{load_data}}
\code{\link{euclidean.clust}}
\code{\link{hamming.clust}}
\code{\link{pearson.clust}}
\code{\link{densitycut.clust}}
\code{\link{evaluate.epiclomal}}
\code{\link{hamming.dist}}
}
