% file MASS/rms.curv.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{rms.curv}
\alias{rms.curv}
\title{
Relative Curvature Measures for Non-Linear Regression
}
\description{
Calculates the root mean square parameter effects and intrinsic relative
curvatures, c^theta and c^iota, for a fitted nonlinear regression, as
defined in Bates \& Watts, section 7.3, p. 253 et seq.
}
\usage{
rms.curv(obj, fit.val=<<see below>>, data=obj$call$data)
}
\arguments{
\item{obj}{
Fitted model object of class \code{nls}.  The model must be fitted using the
default algorithm.
}
\item{fit.val}{
An optional fitted values vector with the gradient matrix and Hessian array
as attributes, as produced by the model function obtained by
using the function \code{deriv3} of David Smith.  Extracted from the
fitted model object by default.
}
\item{data}{
Optional data frame for variables.  Extracted from the fitted model object
call (if any data frame is specified) by default.
}}
\value{
A list of class \code{rms.curv} with components \code{pc} and \code{ic} for parameter
effects and intrinsic relative curvatures multiplied by sqrt(F), \code{ct} and
\code{ci} for c^theta and c^iota (unmultiplied), and \code{C} the C-array as used in
section 7.3.1 of Bates \& Watts.
}
\details{
The method of section 7.3.1 of Bates \& Watts is implemented.  The
function \code{deriv3} should be used generate a model function with first
derivative (gradient) matrix and second derivative (Hessian) array
attributes.  This function should then be used to fit the nonlinear
regression model.


A print method, \code{print.rms.curv}, prints the \code{pc} and \code{ic} components
only, suitably annotated.


If either \code{pc} or \code{ic} exceeds some threshold (0.3 has been suggested) the
curvature is unacceptably high for the planar assumption.
}
\references{
Bates, D. M, and Watts, D. G. (1988) 
\emph{Nonlinear Regression Analysis and its Applications.}
Wiley, New York.
}
\examples{
\dontrun{### Not usable in R
> # The treated sample from the Puromycin data 
> mmcurve <- deriv3(~ Vm * conc/(K + conc), c("Vm", "K"), 
+   function(Vm, K, conc) NULL)
> Treated <- Puromycin[Puromycin$state == "treated", ]
> Purfit1 <- nls(vel ~ mmcurve(Vm, K, conc), data=Treated,
+   start=list(Vm=200, K=0.1))
> rms.curv(Purfit1)
Parameter effects: c^theta x sqrt(F) = 0.2121 
        Intrinsic: c^iota  x sqrt(F) = 0.092 
}}
\keyword{nonlinear}
