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
theta.md(y, u, dfr, limit = 20, eps = .Machine$double.eps^0.25)
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
\code{\link{glm.nb}}, \code{\link{theta.ml}}, \code{\link{theta.mm}}
}
\examples{
\dontrun{theta <- theta.md(y, fitted(fm), dfr = 123)}
}
\keyword{models}
