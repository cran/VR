% file spatial/ppgetregion.d copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{ppgetregion}
\alias{ppgetregion}
\title{
Get Domain for Spatial Point Pattern Analyses
}
\description{
Retrieves the rectangular domain \code{(xl, xu)}  x  \code{(yl, yu)} from the
underlying \code{C} code.
}
\usage{
ppgetregion()
}
\value{
A vector of length four with names `c("xl", "xu", "yl", "yu")'.
}
\seealso{
\code{\link{ppregion}}
}
\keyword{spatial}
