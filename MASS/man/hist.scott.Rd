\name{hist.scott}
\alias{hist.scott}
\alias{hist.FD}
\alias{nclass.scott}
\alias{nclass.FD}
\title{
Plot a Histogram with Automatic Bin Width Selection
}
\description{
Plot a histogram with automatic bin width selection, using the Scott
or Freedman-Diaconis formulae.
}
\usage{
hist.scott(x, prob = TRUE, \dots)
hist.FD(x, prob = TRUE,  \dots)
nclass.scott(x)
nclass.FD(x)
}
\synopsis{
hist.scott(x, prob = TRUE, xlab = deparse(substitute(x)), ...)
hist.FD(x, prob = TRUE, xlab = deparse(substitute(x)), ...)
nclass.scott(x)
nclass.FD(x)
}
\arguments{
\item{x}{
A data vector
}
\item{prob}{
Should the plot have unit area, so be a density estimate?
}
\item{...}{
Further arguments to \code{hist}.
}}
\value{
For the \code{nclass.*} functions, the suggested number of classes.
}
\section{Side Effects}{
The \code{hist.*} function plot a histrogram.
}
\details{
\code{nclass.scott} uses Scott's choice for a normal distribution based on
the estimate of the standard error, and \code{nclass.FD} uses the
Freedman-Diaconis choice based on the inter-quartile range.
}
\references{
Venables, W. N. and Ripley, B. D. (1999)
\emph{Modern Applied Statistics with S-PLUS.}
Springer, pages 118--9.
}
\seealso{
\code{\link{hist}}
}
\keyword{hplot}
\keyword{dplot}
