% file MASS/vcov.nlregb.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{vcov.nlregb}
\alias{vcov.nlregb}
\title{
Calculate Variance-Covariance Matrix for a `nlregb' Fit
}
\description{
Method for \code{vcov} to find the variance-covariance matrix of the
coefficients of fits by \code{nlregb}.
}
\usage{
vcov.nlregb(object, method=c("Fisher", "observed", "Huber"),  
   scale=object$scale, eps=0.001, tol=1)
}
\arguments{
\item{object}{
The return values of a \code{nlregb} fit.
}
\item{method}{
The theoretical basis for the estimate.  This can be based on the
Fisher information (the usual assumption) or the observed information
assuming the model is true or a Huber-White \code{sandwich} estimator
which allows the model to be false. Only the Fisher method is
available if the model was fitted without gradient (Jacobian) information.
}
\item{scale}{
An initial scaling for the parameters.
}
\item{eps}{
The step size (as a multiple of \code{min(1, abs(param))}) for
finite-difference approximations to terms in the Hessian.
}
\item{tol}{
Relative change in sum-of-squares sought in a local quadratic
approximation. See the code for the scaling used.
}}
\value{
A matrix of the estimated covariances between the parameter estimates
in the non-linear regression.
}
\seealso{
\code{\link{vcov}}
}
\keyword{nonlinear}
