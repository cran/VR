% file MASS/orlm.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{orlm}
\alias{orlm}
\alias{orlm.default}
\alias{orlm.formula}
\alias{print.orlm}
\alias{summary.orlm}
\alias{wt.huber}
\alias{hsreg}
\title{
Fit Robust Linear Regression Model
}
\description{
Fits a robust linear regression model, using an M-estimator with
Huber's psi function.
}
\usage{
orlm(formula, data, weights, subset, na.action, 
     model=FALSE, k=1.345, sw=1000, \dots)
}
\arguments{
\item{formula}{
a formula object, with the response on the left of a \code{~} operator,
and the terms, separated by \code{+} operators, on the right.
}
\item{data}{
an optional \code{data.frame} in which to interpret the variables named in the
formula, or in the \code{subset} and the \code{weights} argument.
}
\item{weights}{
optional weights; if supplied, the algorithm fits to minimize the sum
of the weights multiplied into the squared residuals.
The weights must be strictly positive.
}
\item{subset}{
optional expression saying that only a subset of the rows of the data
should be used in the fit.
}
\item{na.action}{
a missing-data filter function, applied to the \code{model.frame}, after
any subset argument has been used.
}
\item{model}{
flag to control what is returned.  If this is \code{TRUE}, then the model frame
is returned.  \code{X} and \code{y} are always returned.
}
\item{k}{
The control value for Winsorizing. The default gives 95\% efficiency at the
normal.
}
\item{sw}{
switch to Huber proposal 2 scale at iteration \code{sw} and beyond.
}
\item{...}{
additional arguments for the fitting routines.
The most likely one is \code{maxit}, which sets the iteration limit, by default
20.
}}
\value{
an object of class \code{"orlm"} representing the fit, inheriting from \code{lm}.
This has all the components of an \code{lm} object, plus \code{k,} the scale \code{s} and
\code{conv} which is a vector monitoring the convergence.
}
\details{
The fit uses Huber's M-estimator, and initially uses the median
absolute deviation scale estimate based on the residuals. This can be
changed to Huber's proposal 2 after \code{sw} iterations.


Generic functions such as \code{print} and \code{summary} have methods to
show the results of the fit.
}
\seealso{
\code{\link{rlm}}
}
\examples{
data(phones)
attach(phones)
res <- orlm(calls ~ year)
print(res)

data(stackloss)
rlm(stack.loss ~ stack.x)
}
\keyword{robust}
