% file MASS/theta.md.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{theta.md}
\alias{theta.md}
\title{
Estimate theta of the Negative Binomial by Deviance
}
\description{
Given the estimated mean vector, estimate \code{theta} of the Negative Binomial
Distribution by equating the deviance to the residual degrees of freedom.
Analogue of a moment estimator.
}
\usage{
theta.md(y, u, dfr, limit=20, eps=sqrt(.Machine$single.eps))
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
\code{\link{glm.nb}}
}
\examples{
\dontrun{theta <- theta.md(y, fitted(fm), dfr=123)}
}
\keyword{models}
% Converted by Sd2Rd version 0.3-2.
