% file MASS/huber.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{huber}
\alias{huber}
\title{
Huber M-estimator of Location with MAD Scale
}
\description{
Finds the Huber M-estimator of location with MAD scale.
}
\usage{
huber(y, k = 1.5, tol = 1e-06)
}
\arguments{
\item{y}{
vector of data values
}
\item{k}{
Winsorizes at \code{k} standard deviations
}
\item{tol}{
convergence tolerance
}}
\value{
list of location and scale parameters

\item{mu}{
location estimate
}
\item{s}{
MAD scale estimate
}}
\references{
Huber, P. J. (1981)
\emph{Robust Statistics.}
Wiley.
}
\seealso{
\code{\link{hubers}}, \code{\link{mad}}
}
\examples{
data(chem)
huber(chem)
}
\keyword{robust}
