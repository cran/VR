\name{glmmPQL}
\alias{glmmPQL}
\title{Fit Generalized Linear Mixed Models via PQL}
\description{
  Fit a GLMM model with multivariate normal random effects, using PQL.
}
\usage{
glmmPQL(fixed, random, family, data, correlation, weights,
        control, niter = 10, verbose = TRUE, ...)
}
\arguments{
  \item{fixed}{a two-sided linear formula object describing the
    fixed-effects part of the model.}
  \item{random}{A formula or list of formulae describing the random effects.}
  \item{family}{a GLM family.}
  \item{data}{an optional data frame used as the first place to find
    variables in the formulae.}
  \item{correlation}{an optional correlation structure.}
  \item{weights}{optional case weights as in \code{\link{glm}}.}
  \item{control}{an optional argument to be passed to
    \code{\link[nlme]{lme}}. (This is modified in the S version by
    passed unchanged in the \R version.)}
  \item{niter}{maximum number of iterations.}
  \item{verbose}{logical: print out record of iterations?}
  \item{\dots}{Further arguments for \code{\link[nlme]{lme}}.}
}
\details{
  \code{glmmPQL} is a wrapper to iterative calls to
  \code{\link[nlme]{lme}}.  Initial values are found by a call to
  \code{glm}, and then \code{lme} is applied to the working vector,
  including the BLUPs of the random effects in forming the linear
  predictor.
}
\value{
  A object of class \code{"lme"}: see \code{\link[nlme]{lmeObject}}.
}
\seealso{\code{\link[nlme]{lme}}}
\references{
  Schall, R. (1991) Estimation in generalized linear models with
  random effects. \emph{Biometrika} \bold{78}, 719--727.

  Breslow, N. E. and Clayton, D. G. (1993) Approximate inference in
  generalized linear mixed models.
  \emph{Journal of the American Statistical Association} \bold{88},
  9--25.

  Wolfinger, R. and O'Connell, M. (1993) Generalized linear mixed models: a
  pseudo-likelihood approach.
  \emph{Journal of Statistical Computation and Simulation} \bold{48}, 233--243.
}
\author{B. D. Ripley}
\examples{
require(nlme)

data(bacteria)
summary(glmmPQL(y ~ trt + I(week> 2), random = ~ 1 | ID, 
                  family = binomial, data = bacteria))
}
\keyword{models}
