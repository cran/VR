% file spatial/pplik.d copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{pplik}
\alias{pplik}
\title{
Pseudo-likelihood Estimation of a Strauss Spatial Point Process
}
\description{
Pseudo-likelihood estimation of a Strauss spatial point process.
}
\usage{
pplik(pp, R, ng=25)
}
\arguments{
\item{pp}{
a pp object
}
\item{R}{
the fixed parameter \code{R}
}
\item{ng}{
use a \code{ng}  x  \code{ng} grid with border \code{R} in the domain for numerical integration.
}}
\value{
estimate for \code{c} in the interval \code{[0, 1]}.
}
\seealso{
\code{\link{Strauss}}
}
\examples{
pines <- ppinit("pines.dat")
pplik(pines, 0.7)
}
\keyword{spatial}
