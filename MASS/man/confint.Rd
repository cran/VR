% file MASS/confint.d
% copyright (C) 1999 W. N. Venables and B. D. Ripley
%
\name{confint}
\alias{confint}
\title{
Confidence Intervals for Model Parameters
}
\description{
Computes confidence intervals for one or more parameters in a fitted
model. Methods currently exist for \code{glm} and \code{nls} fits.
}
\usage{
confint(object, parm = <<see below>>, level = 0.95, \dots)
}
\arguments{
\item{object}{
a fitted model object. Methods currently exist for the classes
\code{"glm"}, \code{"nls"} and for profile objects from these classes.
}
\item{parm}{
a specification of which parameters are to be given confidence
intervals, either a vector of numbers or a vector of names. If
missing, all parameters are considered.
}
\item{level}{
the confidence level required.
}
\item{...}{
additional argument(s) for methods -- the \code{glm} methods allows \code{trace}.
}}
\value{
A matrix (or vector) with columns giving lower and upper confidence
limits for each parameter. These will be labelled as (1-level)/2 and
1 - (1-level)/2 in \% (by default 2.5\% and 97.5\%).
}
\details{
\code{confint} calls the appropriate profile method, then finds the
confidence intervals by interpolation in the profile traces.
If the profile object is already available it should be used as the 
main argument rather than the fitted model object itself.
}
\seealso{
\code{\link{profile}}
}
\examples{
\dontrun{
expn2 <- deriv(~b0 + b1*((w0 - b0)/b1)^(x/d0), 
         c("b0","b1","d0"), function(b0, b1, d0, x, w0) {})


wtloss.init <- function(obj, w0) {
  p <- coef(obj)
  d0 <-  - log((w0 - p["b0"])/p["b1"])/log(2) * p["th"]
  c(p[c("b0", "b1")], d0 = as.vector(d0))
}


out <- NULL
w0s <- c(110, 100, 90)
for(w0 in w0s) {
    fm <- nls(Weight ~ expn2(b0, b1, d0, Days, w0),
              wtloss, start = wtloss.init(wtloss.gr, w0))
    out <- rbind(out, c(coef(fm)["d0"], confint(fm, "d0")))
  }
dimnames(out) <- list(paste(w0s, "kg:"),  c("d0", "low", "high"))
out
}}
\keyword{models}
