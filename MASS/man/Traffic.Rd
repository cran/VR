% file MASS/Traffic.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{Traffic}
\alias{Traffic}
\title{
Effect of Swedish Speed Limits on Accidents
}
\description{
An experiment was performed in Sweden in 1961-2 to assess the
effect of a speed limit on the motorway accident rate.  The
experiment was conducted on 92 days in each year, matched so that
day \code{j} in 1962 was comparable to day \code{j} in 1961.  On some days
the speed limit was in effect and enforced, while on other days
there was no speed limit and cars tended to be driven faster.
The speed limit days tended to be in contiguous blocks. 
}
\format{
This data frame contains the following columns:
}
\arguments{
\item{year}{
1961 or 1962
}
\item{day}{
of year
}
\item{limit}{
was there a speed limit?
}
\item{y}{
traffic accident count for that day
}}
\source{
Svensson, A. (1981)
On the goodness-of-fit test for the multiplicative Poisson model.
\emph{Annals of Statistics,}
\bold{9} 697-704.
}
\keyword{datasets}
