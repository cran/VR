% file MASS/digamma.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{digamma}
\alias{digamma}
\title{
Digamma Function
}
\description{
Calculates the digamma function, which is the first derivative of
the logarithm of the gamma function, for complex vector arguments.
}
\usage{
digamma(z)
}
\arguments{
\item{z}{
Vector argument.  May be real or complex but must have positive
real part.
}}
\value{
Vector of digamma(z) values.
}
\details{
Uses an asymptotic expansion where Re(z) > 5 and a recurrence
formula to such a case where Re(z) <= 5.  Accuracy is good.
}
\references{
M. Abramowitz and I. A. Stegun (1965)
\emph{Handbook of Mathematical Functions.}
Dover.
}
\seealso{
\code{\link{trigamma}}
}
\keyword{math}
