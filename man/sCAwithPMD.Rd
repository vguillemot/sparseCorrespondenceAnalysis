% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sCAwithPMD.R
\name{sCAwithPMD}
\alias{sCAwithPMD}
\title{Sparse Correspondence Analysis based on the "projected Penalized Matrix Decomposition"}
\usage{
sCAwithPMD(
  DATA,
  dimensions = 2L,
  doublecentering = TRUE,
  s1 = rep(1, dimensions),
  s2 = rep(1, dimensions)
)
}
\arguments{
\item{DATA}{the I times J contingency table}

\item{dimensions, }{integer, the number of dimensions to return (default to 2)}

\item{doublecentering, }{logical: should the data be double-centered (default to TRUE)}

\item{s1, }{vector of size 'dimensions" containing the left side regularization parameters, the coefficients in this vector should belong to the interval \link{1, sqrt(I)}, (default to rep(1, dimensions))}

\item{s2, }{vector of size 'dimensions" containing the left side regularization parameters, the coefficients in this vector should belong to the interval \link{1, sqrt(J)}, (default to rep(1, dimensions))}
}
\value{
An object containing all the necessary outputs for sparse CA
}
\description{
Sparse Correspondence Analysis based on the "projected Penalized Matrix Decomposition"
}
\examples{
sCAwithPMD(HairEyeColor[,,1],dimensions = 3L, s1 = rep(0.5 * sqrt(dim(HairEyeColor)[1]), 3), rdsRight = rep(0.5 * sqrt(dim(HairEyeColor)[2]), 3))

}
