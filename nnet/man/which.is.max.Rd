% file nnet/which.is.max.d copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{which.is.max}
\alias{which.is.max}
\title{
Find Maximum Position in Vector
}
\description{
Find the maximum position in a vector, breaking ties at random.
}
\usage{
which.is.max(x)
}
\arguments{
\item{x}{
a vector
}}
\value{
index of a maximal value.
}
\details{
Ties are broken at random.
}
\seealso{
\code{\link{max.col}} in library MASS.
}
\examples{
\dontrun{
pred <- predict(nnet, test)
table(true, apply(pred,1,which.is.max))
}}
\keyword{utilities}
