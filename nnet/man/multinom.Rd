% file nnet/multinom.d copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{multinom}
\alias{multinom}
\title{
Fit Multinomial Log-linear Models
}
\description{
Fits multinomial log-linear models via neural networks.
}
\usage{
multinom(formula, data=sys.parent(), weights, subset, na.action, 
contrasts=NULL, Hess=FALSE, summ=0, censored=FALSE, \dots)
}
\arguments{
\item{formula}{
a formula expression as for regression models, of the form 
\code{response ~ predictors}. The response should be a factor or a matrix with
K columns, which will be interpreted as counts for each of K classes.
A log-linear model is fitted, with coefficients zero for the first
class. An offset can be included: it should be a matrix with K columns
if the response is a matrix with K columns or a factor with K > 2
classes, or a vector for a factor with 2 levels.
See the documentation of \code{formula} for other details.
}
\item{data}{
an optional data frame in which to interpret the variables occurring
in \code{formula}.
}
\item{weights}{
optional case weights in fitting.
}
\item{subset}{
expression saying which subset of the rows of the data should  be used
in the fit. All observations are included by default.
}
\item{na.action}{
a function to filter missing data.
}
\item{contrasts}{
a list of contrasts to be used for some or all of
the factors appearing as variables in the model formula.
}
\item{Hess}{
logical for whether the Hessian (the observed information matrix)
should be returned.
}
\item{summ}{
integer; if non-zero summarize by deleting duplicate rows and adjust weights.
Methods 1 and 2 differ in speed (2 uses \code{C}); method 3 also combines rows
with the same X and different Y, which changes the baseline for the
deviance.
}
\item{censored}{
If Y is a matrix with \code{K > 2} columns, interpret the entries as one
for possible classes, zero for impossible classes, rather than as
counts.
}
\item{...}{
additional arguments for \code{nnet}
}}
\value{
A \code{nnet} object with additional structure.

\item{deviance}{
the residual deviance.
}
\item{edf}{
the (effective) number of degrees of freedom used by the model
}
\item{AIC}{
the AIC for this fit.
}
\item{Hessian}{
(if \code{Hess} is true).
}}
\seealso{
\code{\link{nnet}}
}
\examples{
options(contrasts = c("contr.treatment", "contr.poly"))
library(MASS)
example(birthwt)
bwt.mu <- multinom(low ~ ., bwt)
bwt.mu
\dontrun{Call:
multinom(formula = low ~ ., data = bwt)

Coefficients:
 (Intercept)         age         lwt raceblack raceother  
    0.823477 -0.03724311 -0.01565475  1.192371 0.7406606
     smoke      ptd        ht        ui       ftv1     ftv2+ 
  0.7555234 1.343648 1.913213 0.6802007 -0.4363238 0.1789888

Residual Deviance: 195.4755 
AIC: 217.4755 
}}
\keyword{neural}
\keyword{models}
% Converted by Sd2Rd version 0.3-2.
