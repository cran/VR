% file MASS/logtrans.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{logtrans}
\alias{logtrans}
\alias{logtrans.formula}
\alias{logtrans.lm}
\alias{logtrans.default}
\title{
Estimate log Transformation Parameter
}
\description{
Find and optionally plot the marginal likelihood for alpha
for a transformation model of the form \code{log(y + alpha) ~ x1 + x2 + \dots}.
}
\usage{
logtrans(object, \dots, alpha = seq(0.5, 6, by = 0.25) - min(y),
         plotit = <<see below>>, interp = <<see below>>,
         xlab = "alpha", ylab = "log Likelihood")
}
\synopsis{
logtrans(object, ...)
}
\arguments{
\item{object}{
Fitted linear model object, or formula defining the untransformed
model that is \code{y ~ x1 + x2 + \dots}.  The function is generic.
}
\item{...}{
If \code{object} is a formula, this argument may specify a data frame
as for \code{lm}.
}
\item{alpha}{
Set of values for the transformation parameter, alpha.
}
\item{plotit}{
Should plotting be done?  (Default is \code{TRUE} if a non-null device is
currently active, else \code{FALSE}.)
}
\item{interp}{
Should the marginal log-likelihood be interpolated with a spline
approximation?   (Default is \code{TRUE} if plotting is to be done and
the number of real points is less than 100.)
}
\item{xlab}{
as for \code{plot}.
}
\item{ylab}{
as for \code{plot}.
}}
\value{
List with components \code{x} (for alpha) and \code{y} (for the marginal
log-likelihood values).
}
\section{Side Effects}{
A plot of the marginal log-likelihood is produced, if requested,
together with an approximate mle and 95\% confidence interval.
}
\references{
Venables & Ripley, Chapter 6.
}
\seealso{
\code{\link{boxcox}}
}
\examples{
logtrans(Days ~ Age*Sex*Eth*Lrn, data = quine,
         alpha = seq(0.75, 6.5, len=20))
}
\keyword{regression}
\keyword{models}
\keyword{hplot}
