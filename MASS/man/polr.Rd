% file MASS/polr.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{polr}
\alias{polr}
\alias{extractAIC.polr}
\alias{predict.polr}
\alias{print.polr}
\alias{summary.polr}
\alias{print.summary.polr}
\alias{vcov.polr}
\alias{model.frame.polr}

\title{
Proportional Odds Logistic Regression
}
\description{
Fits a proportional odd logistic regression model to an ordered factor
response.
}
\usage{
polr(formula, data=NULL, weights, start, \dots, subset, 
     na.action=na.fail, contrasts=NULL, Hess=FALSE)
}
\arguments{
\item{formula}{
a formula expression as for regression models, of the form 
\code{response ~ predictors}. The response should be a factor (preferably an ordered
factor), which will be interpreted as an ordinal response, with levels
ordered as in the factor.  A proportional odds model will be
fitted. The model must have an intercept: attempts to remove one will
lead to a warning and be ignored. An offset may be used. See the
documentation of \code{formula} for other details.
}
\item{data}{
an optional data frame in which to interpret the variables occurring
in \code{formula}.
}
\item{weights}{
optional case weights in fitting.
}
\item{start}{
initial values for the parameters.
}
\item{...}{
additional arguments (currently ignored).
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
}}
\value{
A object of class \code{"polr"}.

\item{coefficients}{
the coefficients of the linear predictor.
}
\item{zeta}{
the intercepts for the class boundaries.
}
\item{deviance}{
the residual deviance.
}
\item{fitted.values}{
a matrix, with a column for each level of the response.
}
\item{lev}{
the names of the response levels.
}
\item{terms}{
the \code{terms} structure describing the model.
}
\item{df.residual}{
the number of residual degrees of freedoms, calculated using the weights.
}
\item{edf}{
the (effective) number of degrees of freedom used by the model
}
\item{n}{
the (effective) number of observations,  calculated using the weights
}
\item{call}{
the matched call.
}
\item{convergence}{
the convergence message returned by \code{nlminb}.
}
\item{niter}{
the number of iterations, function and gradient evaluations used by \code{nlminb}.
}
\item{Hessian}{
(if \code{Hess} is true).
}}
\seealso{
\code{\link{nlm}}, \code{\link{glm}}, \code{\link[nnet]{multinom}}.
}
\examples{
options(contrasts=c("contr.treatment", "contr.poly"))
data(housing)
house.plr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
house.plr
summary(house.plr)
predict(house.plr,housing, type="p")
addterm(house.plr, ~.^2, test = "Chisq")
house.plr2 <- stepAIC(house.plr, ~.^2)
house.plr2$anova
}
\keyword{model}
