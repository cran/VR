% file MASS/stepAIC.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{stepAIC}
\alias{stepAIC}
\title{
Choose a model by AIC in a Stepwise Algorithm
}
\description{
Performs stepwise model selection by exact AIC.
}
\usage{
stepAIC(object, scope, scale, direction=c("both", "backward", "forward"), 
        trace=1, keep=NULL, steps=1000, use.start=FALSE, k=2, \dots)
}
\arguments{
\item{object fit}{
an object representing a model of an appropriate class.
This is used as the initial model in the stepwise search.
}
\item{scope}{
defines the range of models examined in the stepwise search. 
}
\item{scale}{
used in the definition of the AIC statistic for selecting the models,
currently only for \code{lm}, \code{aov} and \code{glm} models.
}
\item{direction}{
the mode of stepwise search, can be one of \code{"both"}, \code{"backward"}, or
\code{"forward"}, with a default of \code{"both"}.  If the \code{scope} argument is
missing, the default for \code{direction} is \code{"backward"}.
}
\item{trace}{
if positive, information is printed during the running of \code{stepAIC()}.
Larger values may give more information on the fitting process.
}
\item{keep}{
a filter function whose input is a fitted model object and the 
associated \code{AIC} statistic, and whose output is arbitrary. 
Typically \code{keep} will select a subset of the components of 
the object and return them. The default is not to keep anything.
}
\item{steps}{
the maximum number of steps to be considered.  The default is 1000
(essentially as many as required).  It is typically used to stop the
process early.
}
\item{use.start}{
if true the updated fits are done starting at the linear predictor for
the currently selected model. This may speed up the iterative
calculations for \code{glm} (and other fits), but it can also slow them
down.
}
\item{k}{
the multiple of the number of degrees of freedom used for the penalty.
Only \code{k=2} gives the genuine AIC: \code{k = log(n)} is sometimes referred
to as BIC or SBC.
}
\item{...}{
any additional arguments to \code{extractAIC}. (None are currently used.)
}}
\value{
the stepwise-selected model is returned, with up to two additional
components.  There is an \code{"anova"} component corresponding to the
steps taken in the search, as well as a \code{"keep"} component if the
\code{keep=} argument was supplied in the call. The \code{"Resid. Dev"} column
of the analysis of deviance table refers to a constant minus twice the
maximized log likelihood: it will be a deviance only in cases where a
saturated model is well-defined (thus excluding \code{lm}, \code{aov} and
\code{survreg} fits, for example).
}
\details{
\code{stepAIC} differs from \code{step} and especially \code{step.glm} in using the
exact AIC rather than potentially misleading one-step approximations.
It is also much more widely applicable: all that is required is a
method for \code{extractAIC}, which should return a vector 
\code{c(modeldf, AIC)}.
The default method handles linear models (\code{lm}, \code{aov} and
\code{glm} of family \code{"Gaussian"} with identity link) using \code{addterm.lm} and
\code{dropterm.lm}: for these the results are similar to \code{step.glm} except
that the AIC quoted is Akaike's not Hastie's. (The additive constant
is chosen so that in that case AIC is identical to Mallows' Cp if the
scale is known.)


There is a potential problem in using \code{glm} fits with a variable
\code{scale}, as in that case the deviance is not simply related to the
maximized log-likelihood. The function \code{extractAIC.glm} makes the
appropriate adjustment for a \code{gaussian} family, but may need to be
amended for other cases. (The \code{binomial} and \code{poisson} families have
fixed \code{scale} by default and do not correspond to a particular
maximum-likelihood problem for variable \code{scale}.)


Where a conventional deviance exists (e.g. for \code{lm}, \code{aov} and \code{glm}
fits) this is quoted in the analysis of variance table: it is the
\emph{unscaled} deviance.
}
\seealso{
\code{\link{addterm}}, \code{\link{dropterm}}, \code{\link{step}}
}
\examples{
data(quine)
quine.hi <- aov(log(Days + 2.5) ~ .^4, quine)
quine.nxt <- update(quine.hi, . ~ . - Eth:Sex:Age:Lrn)
quine.stp <- stepAIC(quine.nxt, 
    scope = list(upper = ~Eth*Sex*Age*Lrn, lower = ~1), 
    trace = FALSE)
quine.stp$anova

data(cpus)
cpus1 <- cpus
attach(cpus)
for(v in names(cpus)[2:7]) 
  cpus1[[v]] <- cut(cpus[[v]], unique(quantile(cpus[[v]])), 
                    include.lowest = TRUE)
detach()
cpus0 <- cpus1[, 2:8]  # excludes names, authors' predictions
cpus.samp <- sample(1:209, 100)
cpus.lm <- lm(log10(perf) ~ ., data=cpus1[cpus.samp,2:8])
cpus.lm2 <- stepAIC(cpus.lm, trace=FALSE)
cpus.lm2$anova

example(birthwt)
birthwt.glm <- glm(low ~ ., family=binomial, data=bwt)
birthwt.step <- stepAIC(birthwt.glm, trace=FALSE)
birthwt.step$anova
birthwt.step2 <- stepAIC(birthwt.glm, ~ .^2 + I(scale(age)^2)
    + I(scale(lwt)^2), trace=FALSE)
birthwt.step2$anova


quine.nb <- glm.nb(Days ~ .^4, data=quine)
quine.nb2 <- stepAIC(quine.nb)
quine.nb2$anova
}
\keyword{models}
