\name{MASS-internal}
\alias{biplot.bdr}
\alias{enlist}
\alias{fbeta}
\alias{frequency.polygon}
\alias{nclass.freq}
\alias{neg.bin}
\alias{negexp.SSival}
\alias{pairs.profile}
\alias{plot.profile}
\alias{print.Anova}
\alias{print.abbrev}
\title{Internal MASS functions}
\description{
 Internal MASS functions
}
\usage{
biplot.bdr(obs, bivars, col, cex = rep(par("cex"), 2),
           olab = NULL, vlab = NULL, xlim = NULL, ylim = NULL, \dots)
enlist(vec)
fbeta(x, alpha, beta)
frequency.polygon(x, nclass = nclass.freq(x), xlab="", ylab="", \dots)
nclass.freq(x)
neg.bin(theta = stop("theta must be given"))
negexp.SSival(mCall, data, LHS)
pairs.profile(x, colours = 2:3, \dots)
plot.profile(x, nseg, \dots)
print.Anova(x, \dots)
}
\details{
  These are not to be called by the user.  Some are for compatibility
  with earlier versions of MASS (the book).
}
\keyword{internal}
