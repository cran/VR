% file MASS/isoMDS.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{isoMDS}
\alias{isoMDS}
\alias{Shepard}
\title{
Kruskal's Non-metric Multidimensional Scaling
}
\description{
One form of non-metric multidimensional scaling
}
\usage{
isoMDS(d, y = cmdscale(d, k), k = 2, maxit = 50, trace = TRUE)
}
\arguments{
\item{d}{
  distance structure of the form returned by \code{dist}, or a full,
  symmetric matrix.  Data are assumed to be dissimilarities or relative
  distances, but must be positive except for self-distance.
}
\item{y}{
  An initial configuration. If none is supplied, \code{cmdscale} is used
  to provide the classical solution. The dimension of the initial
  configuration determines that of the answer.
}
\item{k}{
  The desired dimension for the solution, passed to \code{cmdscale}.
}
\item{maxit}{
  The maximum number of iterations.
}
\item{trace}{
  Logical for tracing optimization. Default \code{true}.
}}
\value{
Two components:

\item{points}{
  A two-column vector of the fitted configuration.
}
\item{stress}{
  The final stress achieved (in percent).
}}
\section{Side Effects}{
  If \code{trace} is true, the initial stress and the current stress
  are printed out every 5 iterations.
}
\details{
This chooses a k-dimensional (default k = 2) configuration to minimize
the stress, the square root of the ratio of the the sum of squared
differences between the input distances and those of the configuration
to the sum of configuration distances squared. However, the input
distances are allowed a monotonic transformation.

An iterative algorithm is used, which will usually converge in around
10 iterations. As this is necessarily an \eqn{O(n^2)} calculation, it is slow
for large datasets. Further, since the configuration is only determined
up to rotations and reflections (by convention the centroid is at the
origin), the result can vary considerably from machine to machine.
}
\references{
T. F. Cox and M. A. A. Cox (1994, 2001)
\emph{Multidimensional Scaling}.
Chapman & Hall.
}
\seealso{
\code{\link{cmdscale}}, \code{\link{sammon}}
}
\examples{
library(mva) # for dist
data(swiss)
swiss.x <- as.matrix(swiss[, -1])
swiss.dist <- dist(swiss.x)
swiss.mds <- isoMDS(swiss.dist)
plot(swiss.mds$points, type = "n")
text(swiss.mds$points, labels = as.character(1:nrow(swiss.x)))
swiss.sh <- Shepard(swiss.dist, swiss.mds$points)
plot(swiss.sh, pch = ".")
lines(swiss.sh$x, swiss.sh$yf, type = "S")
}
\keyword{multivariate}
