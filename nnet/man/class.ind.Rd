% file nnet/class.ind.d copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{class.ind}
\alias{class.ind}
\title{
Generates Class Indicator Matrix from a Factor
}
\description{
Generates a class indicator function from a given factor.
}
\usage{
class.ind(cl)
}
\arguments{
\item{cl}{
factor or vector of classes for cases.
}}
\value{
a matrix which is zero except for the column corresponding to the class.
}
\examples{
# The function is currently defined as
function(cl)
{
  n <- length(cl)
  cl <- as.factor(cl)
  x <- matrix(0, n, length(levels(cl)))
  x[(1:n) + n * (codes(cl) - 1)] <- 1
  dimnames(x) <- list(row.names(cl), levels(cl))
  x
}
}
\keyword{neural}
\keyword{utilities}
