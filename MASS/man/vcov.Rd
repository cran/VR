% file MASS/vcov.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{vcov}
\alias{vcov}
\alias{vcov.lm}
\alias{vcov.glm}
\alias{vcov.nls}
\alias{vcov}
\title{
Calculate Variance-Covariance Matrix for a Fitted Model Object
}
\description{
Returns the variance-covariance matrix.
This is a generic function.
Functions with names beginning in \code{vcov.} will be
methods for this function.
Classes with methods for this function include:
\code{lm}, \code{glm}, \code{nls}.
}
\usage{
vcov(object, \dots)
}
\arguments{
\item{object}{
a fitted model object.
}
\item{...}{
additional arguments for method functions.
}}
\value{
a matrix of the estimated covariances between the parameter estimates
in the linear or non-linear predictor of the model.
}
\keyword{models}
\keyword{nonlinear}
