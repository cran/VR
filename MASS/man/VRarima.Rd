\name{VRarima}
\alias{VRarima}
\alias{print.VRarima}
\alias{predict.VRarima}
\alias{diagnostics}
\alias{diagnostics.VRarima}
\title{
Wrapper Functions for ARIMA Modelling
}
\description{
Wrapper functions for ARIMA modelling.  These provide interfaces to
\code{\link[ts]{arima0}}, \code{\link[ts]{predict.arima0}}
and \code{\link[ts]{arima0.diag}}.
}
\usage{
VRarima(x, order, seasonal, n.cond, \dots)
\method{predict}{VRarima}(object, n.ahead = 1, se.fit = TRUE, \dots)
\method{diagnostics}{VRarima}(x, \dots)
}
\arguments{
\item{x}{
a time series.
}
\item{order}{
a vector of three components, (p, d, q).
}
\item{seasonal}{
  list with components \code{order} and \code{period}:
  \code{order} is  vector of three components, (p, d, q).
}
\item{n.cond}{
Ignored: for S-PLUS compatibility.
}
\item{object}{
an object of class \code{"VRarima"}.
}
\item{n.ahead}{
The number of steps ahead to predict.
}
\item{se.fit}{
logical: should standard errors be returned?
}
\item{\dots}{
  further arguments to \code{arima0}, \code{predict.arima0} or
  \code{arima0.diag}.
}}
\value{
  For \code{VRarima}, an object of class \code{"VRarima"}
  which inherits from \code{"arima0"}.

  For \code{predict.VRarima}, a vector if \code{se.fit = FALSE}, or a
  list with components \code{pred} and \code{se}.
}
\seealso{
  \code{\link[ts]{arima0}}, \code{\link[ts]{predict.arima0}} and
  \code{\link[ts]{arima0.diag}}.
}
\keyword{ts}
