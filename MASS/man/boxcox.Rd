% file MASS/boxcox.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{boxcox}
\alias{boxcox}
\alias{boxcox.default}
\alias{boxcox.formula}
\alias{boxcox.lm}
\title{
Box-Cox Transformations for Linear Models
}
\description{
Computes and optionally plots profile log-likelihoods for the
parameter of the Box-Cox simple power transformation \code{y^lambda}.
}
\usage{
boxcox(object, lambda, plotit, interp, eps, xlab, ylab, \dots)
}
\synopsis{
boxcox(object, ...)
}
\arguments{
\item{object}{
a formula or fitted model object.  Currently only \code{lm} and \code{aov} objects
are handled.
}
\item{lambda}{
vector of values of \code{lambda} -- default \eqn{(-2, 2)} in steps of 0.1
}
\item{plotit}{
logical which controls whether the result should be plotted.
Default to \code{true} if a graphics device is currently open.
}
\item{interp}{
logical which controls whether spline interpolation is used.
Default to \code{true} if plotting with \code{lambda} of length less than 100.
}
\item{eps}{
Tolerance for \code{lambda = 0}; defaults to 0.02
}
\item{xlab}{
defaults to \code{"lambda"}
}
\item{ylab}{
defaults to \code{"log-Likelihood"}
}
\item{\dots}{
additional parameters to be used in the model fitting.
}}
\value{
A list of the \code{lambda} vector and the computed profile
log-likelihood vector, invisibly if the result is plotted.
}
\section{Side Effects}{
If \code{plotit = TRUE} plots loglik \emph{vs} \code{lambda} and indicates a 95\%
confidence interval about the maximum observed value of \code{lambda}. If
\code{interp = TRUE}, spline interpolation is used to give a smoother plot.
}
\references{
  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.} Fourth edition.  Springer.
}
\examples{
data(trees)
boxcox(Volume ~ log(Height) + log(Girth), data = trees,
       lambda = seq(-0.25, 0.25, length = 10))

boxcox(Days+1 ~ Eth*Sex*Age*Lrn, data = quine,
       lambda = seq(-0.05, 0.45, len = 20))
}
\keyword{regression}
\keyword{models}
\keyword{hplot}
