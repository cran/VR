% file nnet/which.is.max.d copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{which.is.max}
\alias{which.is.max}
\alias{max.col}
\title{
Find Maximum Position in Vector or Matrix
}
\description{
Find the maximum position in a vector or for each row of a matrix, breaking ties at random.
}
\usage{
which.is.max(x)
max.col(A)
}
\arguments{
\item{x}{
a vector
}
\item{A}{
a matrix
}}
\value{
index of a maximal value (for each row in the case of \code{max.col}).
}
\details{
Ties are broken at random.
}
\examples{
\dontrun{
pred <- predict(nnet, test)
table(true, apply(pred,1,which.is.max))
}}
\keyword{neural}
\keyword{utilities}
% Converted by Sd2Rd version 0.3-2.
