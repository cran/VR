\name{fitdistr}
\alias{fitdistr}
\alias{print.fitdistr}
\alias{coef.fitdistr}
\title{
Maximum-likelihood Fitting of Univariate Distributions
}
\description{
Maximum-likelihood fitting of univariate distributions, allowing
parameters to be held fixed if desired.
}
\usage{
fitdistr(x, densfun, start, \dots)
}
\arguments{
\item{x}{
A numeric vector.
}
\item{densfun}{
Either a character string or a function returning a density evaluated
at its first argument.


Distributions \code{"beta"}, \code{"cauchy"}, \code{"chi-squared"}, \code{"exponential"},
\code{"f"}, \code{"gamma"}, \code{"log-normal"}, \code{"lognormal"}, \code{"logistic"},
\code{"negative binomial"}, \code{"normal"}, \code{"t"}, \code{"uniform"} and \code{"weibull"}
are recognised, case being ignored.
}
\item{start}{
A named list giving the parameters to be optimized with initial
values.  This can be omitted for some of the named distributions
(see Details).
}
\item{\dots}{
Additional parameters, either for \code{densfun} or for \code{optim}.  In
particular, it can be used to specify bounds via \code{lower} or \code{upper} or
both.  If arguments of \code{densfun} (or the density function
corresponding to a character-string specification) are included they
will be held fixed.
}}
\value{
an object of class \code{"fitdistr"}, a list with two components,

\item{estimate}{
the parameter estimates, and
}
\item{sd}{
the estimated standard errors.
}}
\details{
Direct optimization of the log-likelihood is performed, with numerical
derivatives.  The estimated standard errors are taken from the
observed information, calculated by a numerical approximation.

For the following named distributions, reasonable starting values will
be computed if \code{start} is omitted or only partially specified:
\code{cauchy}, \code{gamma}, \code{logistic}, \code{negative binomial}
( parametrized by
\code{mu} and \code{size}), \code{t}, \code{uniform}, \code{weibull}.
}
\references{
  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.} Fourth edition.  Springer.
}
\examples{
set.seed(123)
x <- rgamma(100, shape = 5, rate = 0.1)
fitdistr(x, "gamma")
## now do this directly with more control.
fitdistr(x, dgamma, list(shape = 1, rate = 0.1), lower = 0.01)


set.seed(123)
x2 <- rt(250, df = 9)
fitdistr(x2, "t", df = 9)
## allow df to vary: not a very good idea!
fitdistr(x2, "t")
## now do this directly with more control.
mydt <- function(x, m, s, df) dt((x-m)/s, df)/s
fitdistr(x2, mydt, list(m = 0, s = 1), df = 9, lower = c(-Inf, 0))


set.seed(123)
x3 <- rweibull(100, shape = 4, scale = 100)
fitdistr(x3, "weibull")


set.seed(123)
x4 <- rnegbin(500, mu = 5, theta = 4)
fitdistr(x4, "Negative Binomial") # R only
}
\keyword{distribution}
\keyword{htest}
