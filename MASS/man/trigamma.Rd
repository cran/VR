% file MASS/trigamma.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{trigamma}
\alias{trigamma}
\title{
Trigamma Function
}
\description{
Calculates the trigamma function, which is the second derivative of
the logarithm of the gamma function, for complex vector arguments.
}
\usage{
trigamma(z)
}
\arguments{
\item{z}{
Vector argument.  May be real or complex but must have positive
real part.
}}
\value{
Vector of trigamma(z) values.
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
\code{\link{digamma}}
}
\keyword{math}
