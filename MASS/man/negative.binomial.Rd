% file MASS/negative.binomial.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{negative.binomial}
\alias{negative.binomial}
\title{
Family function for Negative Binomial GLMs 
}
\description{
Specifies the information required to fit a Negative Binomial generalized
linear model, with known \code{theta} parameter, using \code{glm()}.
}
\usage{
negative.binomial(theta=stop("theta must be specified"), link="log")
}
\arguments{
\item{theta}{
The known value of the additional parameter, \code{theta}.
}
\item{link}{
The link function.  Currently must be one of \code{log}, \code{sqrt} or \code{identity}.
}}
\value{
A list of functions and expressions needed by \code{glm()} to fit a Negative
Binomial generalized linear model.
}
\seealso{
\code{\link{glm.nb}}, \code{\link{anova.negbin}}, \code{\link{summary.negbin}}
}
\examples{
# Fitting a Negative Binomial model to the quine data 
#   with theta=2 assumed known. 
#
data(quine)
quine.nb <- glm(Days ~ .^4, family=neg.bin(2), data=quine)
}
\keyword{regression}
\keyword{models}
