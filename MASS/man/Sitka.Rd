% file MASS/Sitka.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{Sitka}
\alias{Sitka}
\title{
Growth Curves for Sitka Spruce Trees in 1988
}
\description{
The \code{Sitka} data frame has 395 rows and 4 columns.  It gives repeated
measurements on the log-size of 79 Sitka spruce trees, 54 of which
were grown in ozone-enriched chambers and 25 were controls.  The size
was measured five times in 1988, at roughly monthly intervals.
}
\format{
This data frame contains the following columns:
}
\arguments{
\item{size}{
measured size (height times diameter squared) of tree, on log scale.
}
\item{Time}{
time of measurement in  days since 1 January 1988
}
\item{tree}{
number of tree
}
\item{treat}{
either \code{"ozone"} for an ozone-enriched chamber or \code{"control"}
}}
\seealso{
\code{\link{Sitka89}}
}
\source{
P. J. Diggle, K.-Y. Liang and S. L. Zeger (1994)
\emph{Analysis of Longitudinal Data.}
Clarendon Press, Oxford
}
\keyword{datasets}
