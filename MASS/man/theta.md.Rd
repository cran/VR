% file MASS/theta.md.d
% copyright (C) 1994-2004 W. N. Venables and B. D. Ripley
%
\name{theta.md}
\alias{theta.md}
\alias{theta.ml}
\alias{theta.mm}
\title{
Estimate theta of the Negative Binomial
}
\description{
  Given the estimated mean vector, estimate \code{theta} of the
  Negative Binomial Distribution.
}
\usage{
theta.md(y, mu, dfr, limit = 20, eps = .Machine$double.eps^0.25)

theta.ml(y, mu, n, limit = 10, eps = .Machine$double.eps^0.25,
         trace = FALSE)

theta.mm(y, mu, dfr, limit = 10, eps = .Machine$double.eps^0.25)
}
\arguments{
  \item{y}{
    Vector of observed values from the Negative Binomial.
  }
  \item{mu}{
    Estimated mean vector.
  }
  \item{n}{
    Number of data points (defaults to length of \code{y})
  }
  \item{dfr}{
    Residual degrees of freedom (assuming \code{theta} known).
  }
  \item{limit}{
    Limit on the number of iterations.
  }
  \item{eps}{
    Tolerance to determine convergence.
  }
  \item{trace}{
    logical: should iteration progress be printed?
  }
}
\details{
  \code{theta.md} estimates by equating the deviance to the residual
  degrees of freedom, an analogue of a moment estimator.

  \code{theta.ml} uses maximum likelihood.

  \code{theta.mm} calculates the moment estimator of \code{theta} by
  equating the Pearson chi-square
  \eqn{sum((y-\mu)^2/(\mu+\mu^2/\theta))} to the residual
  degrees of freedom.
}
\value{
  The required estimate of \code{theta}, as a scalar.
  For \code{theta.ml}, the standard error is given as attribute \code{"SE"}.
}
\seealso{
  \code{\link{glm.nb}}
}
\examples{
quine.nb <- glm.nb(Days ~ .^2, data = quine)
theta.md(quine$Days, fitted(quine.nb), dfr = df.residual(quine.nb))
theta.ml(quine$Days, fitted(quine.nb))
theta.mm(quine$Days, fitted(quine.nb), dfr = df.residual(quine.nb))
}
\keyword{models}
