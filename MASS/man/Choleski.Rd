% file MASS/Choleski.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{Choleski}
\alias{Choleski}
\title{
Perform Choleski Decomposition of a Positive Definite Symmetric Matrix
}
\description{
Perform Choleski decomposition of a real positive-definite symmetric
Matrix (that was generated using \code{library(Matrix)}).
}
\usage{
Choleski(x)
}
\arguments{
\item{x}{
A real symmetric matrix
}}
\value{
A lower triangular Matrix \code{L} such that \code{L \%*\% t(L) = A}.
}
\details{
An error will result unless \code{A} is strictly positive definite.
}
\seealso{
\code{\link{chol}}
}
\examples{
### Not usable in R
\dontrun{
> library(Matrix)
> A <- Matrix(c(6,3,3,4),2,2)
> Choleski(A)
       [,1]   [,2] 
[1,] 2.4495 0.0000
[2,] 1.2247 1.5811
attr(, "class"):
[1] "LowerTriangular" "Matrix"
}}
\keyword{algebra}
