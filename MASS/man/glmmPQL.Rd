\name{glmmPQL}
\alias{glmmPQL}
\title{
Fit Generalized Linear Mixed Models via PQL
}
\description{
Fit a GLMM model with multivariate normal random effects, using
Penalized Quasi-Likelihood.
}
\usage{
glmmPQL(fixed, random, family, data, correlation, weights,
        control, niter = 10, verbose = TRUE, \dots)
}
\arguments{
\item{fixed}{
a two-sided linear formula giving fixed-effects part of the model.
}
\item{random}{
A formula or list of formulae describing the random effects.
}
\item{family}{
a GLM family.
}
\item{data}{
an optional data frame used as the first place to find
variables in the formulae.
}
\item{correlation}{
an optional correlation structure.
}
\item{weights}{
optional case weights as in \code{glm}.
}
\item{control}{
an optional argument to be passed to \code{lme}.
}
\item{niter}{
maximum number of iterations.
}
\item{verbose}{
logical: print out record of iterations?
}
\item{...}{
Further arguments for \code{lme}.
}}
\value{
A object of class \code{"lme"}: see \code{\link[nlme]{lmeObject}}.
}
\references{
Schall, R. (1991) Estimation in generalized linear models with
random effects.
\emph{Biometrika}
\bold{78}, 719--727.


Breslow, N. E. and Clayton, D. G. (1993) Approximate inference in
generalized linear mixed models.
\emph{Journal of the American Statistical Association}
\bold{88}, 9--25.


Wolfinger, R. and O'Connell, M. (1993) Generalized linear mixed models: a
pseudo-likelihood approach.
\emph{Journal of Statistical Computation and Simulation}
\bold{48}, 233--243.
}
\seealso{
  \code{\link[nlme]{lme}}
}
\examples{
library(nlme)
data(bacteria)
summary(glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
                family = binomial, data = bacteria))
}
\keyword{models}
