% file MASS/theta.mm.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{theta.mm}
\alias{theta.mm}
\title{
Estimate theta of the Negative Binomial by Moments
}
\description{
Given the mean estimate, calculate the moment estimator of theta by
equating
\code{sum((y-mu)^2/(mu+mu^2/theta))} 
to the residual degrees of
freedom.
}
\usage{
theta.mm(y, u, dfr, limit=10, eps=.Machine$double.eps^0.25)
}
\arguments{
\item{y}{
Vector of observed values from the Negative Binomial.
}
\item{u}{
Estimated mean vector.
}
\item{dfr}{
Residual degrees of freedom (assuming \code{theta} known).
}
\item{limit}{
Limit on the number of iterations.
}
\item{eps}{
Tolerance to determine convergence.
}}
\value{
The required estimate of \code{theta}, as a scalar.
}
\seealso{
\code{\link{glm.nb}}, \code{\link{theta.ml}}, \code{\link{theta.md}}
}
\examples{
\dontrun{theta <- theta.mm(y, fitted(fm), dfr=123)}
}
\keyword{models}
