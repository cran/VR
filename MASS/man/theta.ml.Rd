% file MASS/theta.ml.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{theta.ml}
\alias{theta.ml}
\title{
Estimate theta of the Negative Binomial by Maximum Likelihood
}
\description{
Given the estimated mean vector, estimate \code{theta} of the Negative Binomial
Distribution.
}
\usage{
theta.ml(y, mu, n, limit = 10, eps = .Machine$double.eps^0.25, trace = FALSE)
}
\arguments{
\item{y}{
Vector of observed values from the Negative Binomial.
}
\item{mu}{
Estimated mean vector.
}
\item{n}{
Number of data points (defaults to length of \code{y})
}
\item{limit}{
Limit on the number of iterations.
}
\item{eps}{
Tolerance to determine convergence.
}
\item{trace}{
logical: should iteration progress be printed?
}}
\value{
The required estimate of \code{theta}, as a scalar.
}
\seealso{
\code{\link{glm.nb}}, \code{\link{theta.md}}, \code{\link{theta.mm}}
}
\examples{
\dontrun{theta <- theta.ml(y, fitted(fm), dfr = 123)}
}
\keyword{models}
