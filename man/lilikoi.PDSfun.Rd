% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lilikoi.PDSfun.r
\name{lilikoi.PDSfun}
\alias{lilikoi.PDSfun}
\title{A PDSfun Function}
\usage{
lilikoi.PDSfun(qvec)
}
\arguments{
\item{qvec}{This is the Metabolite_pathway_table from MetaTOpathway function. This table includes the metabolites ids and the its corssponding hmdb ids}
}
\value{
A large matrix of the pathway deregulation scores for each pathway in different samples.
}
\description{
This function allows you to compute Pathway Desregulation Score deriving
make sure that you have the below database for the metabolites and pathway list:
meta_path.RData
}
\examples{
\donttest{
dt <- lilikoi.Loaddata(file=system.file("extdata",
  "plasma_breast_cancer.csv", package = "lilikoi"))
Metadata <- dt$Metadata
dataSet <- dt$dataSet
convertResults=lilikoi.MetaTOpathway('name')
Metabolite_pathway_table = convertResults$table
# PDSmatrix= lilikoi.PDSfun(Metabolite_pathway_table)
}
}
\references{
Nygård, S., Lingjærde, O.C., Caldas, C. et al.
PathTracer: High-sensitivity detection of differential pathway activity in tumours.
Sci Rep 9, 16332 (2019). https://doi.org/10.1038/s41598-019-52529-3
}
\keyword{PDS}
