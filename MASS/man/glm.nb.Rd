% file MASS/glm.nb.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{glm.nb}
\alias{glm.nb}
\alias{family.negbin}
\title{
Fit a Negative Binomial Generalized Linear Model
}
\description{
A modification of the system function \code{glm()} to include estimation of the
additional parameter, \code{theta}, for a Negative Binomial generalized linear
model.
}
\usage{
glm.nb(formula, \dots,  init.theta, link=log)
}
\synopsis{
glm.nb(formula = formula(data), data = sys.parent(), weights,
                 subset, na.action, start = eta,
                 control = glm.control(...), method = "glm.fit",
                 model = F, x = F, y = T, contrasts = NULL, ...,
                 init.theta, link = log)
}
\arguments{
\item{formula}{
Model formula, as for \code{glm()}.
}
\item{...}{
Any other arguments for the \code{glm()} function except \code{family}.
}
\item{init.theta}{
Optional initial value for the theta parameter.  If omitted a moment
estimator after an initial fit using a Poisson GLM is used.
}
\item{link}{
The link function.  Currently must be one of \code{log}, \code{sqrt} or \code{identity}.
}}
\value{
A fitted model object of class \code{negbin} inheriting from \code{glm} and \code{lm}.
The object is like the output of \code{glm} but contains three additional
components, namely \code{theta} for the ML estimate of theta, \code{SE.theta} for
its approximate standard error (using observed rather than expected
information), and \code{twologlik} for twice the log-likelihood function.
}
\details{
An alternating iteration process is used.  For given \code{theta} the GLM
is fitted using the same process as used by \code{glm()}.  For fixed means
the \code{theta} parameter is estimated using score and information
iterations.  The two are alternated until convergence of both. (The
number of alternations and the number of iterations when estimating
\code{theta} are controlled by the \code{maxit} parameter of \code{glm.control}.)


Setting \code{trace > 0} traces the alternating iteration process. Setting
\code{trace > 1} traces the \code{glm} fit, and setting \code{trace > 2} traces the
estimation of \code{theta}.
}
\seealso{
\code{\link{glm}}, \code{\link{negative.binomial}}, \code{\link{anova.negbin}}, \code{\link{summary.negbin}},
\code{\link{theta.ml}}, \code{\link{theta.md}}, \code{\link{theta.mm}}
}
\examples{
data(quine)
quine.nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data=quine)
quine.nb2 <- update(quine.nb1, . ~ . + Sex:Age:Lrn)
quine.nb3 <- update(quine.nb2, Days ~ .^4)
anova(quine.nb1, quine.nb2, quine.nb3)
}
\keyword{regression}
\keyword{models}
