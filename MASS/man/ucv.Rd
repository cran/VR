% file MASS/ucv.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{ucv}
\alias{ucv}
\title{
Unbiased Cross-Validation for Bandwidth Selection
}
\description{
Uses unbiased cross-validation to select the bandwidth of a Gaussian
kernel density estimator.
}
\usage{
ucv(x, nb = 1000, lower, upper)
}
\arguments{
\item{x}{
a numeric vector
}
\item{nb}{
number of bins to use.
}
\item{lower, upper}{
Range over which to minimize.  The default is almost always satisfactory.
}}
\value{
a bandwidth.
}
\references{
Scott, D. W. (1992)
\emph{Multivariate Density Estimation: Theory, Practice, and  Visualization.}
Wiley.
}
\seealso{
\code{\link{bcv}}, \code{\link{width.SJ}}, \code{\link{density}}
}
\examples{
data(geyser)
ucv(geyser$duration)
}
\keyword{dplot}
