% file MASS/histplot.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{histplot}
\alias{histplot}
\title{
Trellis of Histograms
}
\description{
Plot a trellis of histograms of data.
}
\usage{
histplot(formula, data, 
         prepanel = prepanel.histplot, panel = panel.histplot, 
         nbins = 5, h, x0 = -h/1000, breaks, prob = TRUE, 
         aspect = "fill", \dots)
}
\arguments{
\item{formula}{
a formula in the form:
~ x | g1 * g2 * \dots{}
where the given variables \code{g1}, \code{g2}, \dots{} may be omitted.
The \code{x} variable is a numeric response and is the data used to
construct each histogram.
Variables on the right hand side of the vertical bar are factors or shingles,
and are used to parcel the \code{x} values out to various panels.
If there are no given variables the initial ~ may be omitted.
}
\item{data}{
A data frame as list against which to evaluate the formula.
}
\item{prepanel}{
The prepanel function to be used.
}
\item{panel}{
The panel function to be used.
}
\item{nbins}{
The suggested number of bins. The default is small, as the number
of data points in each panel is unknown at this stage.
}
\item{h}{
The width of the bins.
}
\item{x0}{
A shift for the bin sequence - the breakpoints are at \code{x0 + n * h}.
}
\item{breaks}{
The set of breakpoints to be used. (Usually omitted).
}
\item{prob}{
If \code{true} (the default) plot a true histogram.
The vertical axis has a 
\emph{relative frequency density}
scale, so the product of the dimensions of any panel gives the
relative frequency.  Hence the total area under the histogram
is 1 and it is directly comparable with most other estimates
of the probability density function.
If \code{false} plot the counts in the bins.
}
\item{...}{
generic Trellis arguments. See under \code{trellis.args} for details.
}}
\value{
an object of class \code{trellis},
which is automatically plotted by \code{print.trellis}.
}
\details{
If \code{breaks} is specified, those breakpoints are used. Otherwise if \code{h}
is specified, a regular grid of bins is used with width \code{h}.
If neither \code{breaks} nor \code{h} is specified, \code{nbins} is used to select
a suitable \code{h}.  Beware of using bandwidth selection algorithms to
choose \code{h}, as the data will be split, possibly unevenly, across panels.


Attempts to use unequally spaced breaks with \code{prob=FALSE} will be honoured
with a warning, but this is strongly discouraged as misleading.
}
\seealso{
\code{\link{truehist}}, \code{\link{histogram}}, \code{\link{trellis.args}}
}
\keyword{trellis}
\keyword{hplot}
\keyword{dplot}
% Converted by Sd2Rd version 0.3-2.
