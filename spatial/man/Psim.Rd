% file spatial/Psim.d copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{Psim}
\alias{Psim}
\title{
Simulate Binomial Spatial Point Process
}
\description{
Simulate Binomial spatial point process.
}
\usage{
Psim(n)
}
\arguments{
\item{n}{
number of points
}}
\value{
list of vectors of \code{x} and \code{y} coordinates.
}
\details{
relies on the region being set by \code{ppinit} or \code{ppregion}.
}
\section{Side Effects}{
uses the random number generator.
}
\seealso{
\code{\link{SSI}}, \code{\link{Strauss}}
}
\examples{
towns <- ppinit("towns.dat")
par(pty="s")
plot(Kfn(towns, 10), type="s", xlab="distance", ylab="L(t)")
for(i in 1:10) lines(Kfn(Psim(69), 10))
}
\keyword{spatial}
