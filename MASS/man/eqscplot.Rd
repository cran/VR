% file MASS/eqscplot.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{eqscplot}
\alias{eqscplot}
\title{
Plots with Geometrically Equal Scales
}
\description{
Version of a scatterplot with scales chosen to be equal on both axes, that
is 1cm represents the same units on each
}
\usage{
eqscplot(x, y, ratio = 1, tol = 0.04, uin, \dots)
}
\synopsis{
eqscplot(x, y, ratio = 1, tol = 0.04, uin, xlim = range(x, na.rm=TRUE),
         ylim = range(y, na.rm=TRUE), xlab, ylab, ...)
}
\arguments{
\item{x}{
vector of x values, or a 2-column matrix, or a list with components
\code{x} and \code{y}
}
\item{y}{
vector of y values
}
\item{ratio}{
desired ratio of units on the axes. Units on the y axis are drawn at
\code{ratio} times the size of units on the x axis. Ignored if \code{uin} is
specified and of length 2.
}
\item{tol}{
proportion of white space at the margins of plot
}
\item{uin}{
desired values for the units per inch parameter. If of length 1, the
desired units per inch on the x axis.
}
\item{...}{
further arguments for \code{plot}
}}
\value{
invisibly, the values of \code{uin} used for the plot.
}
\section{Side Effects}{
performs the plot.
}
\details{
Limits for the x and y axes are chosen so that they include the
data. One of the sets of limits is then stretched from the midpoint to
make the units in the ratio given by \code{ratio}. Finally both are
stretched by \code{1 + tol} to move points away from the axes, and the
points plotted.
}
\section{Note}{
Arguments \code{ratio} and \code{uin} were suggested by Bill Dunlap.
}
\seealso{
\code{\link{plot}}, \code{\link{par}}
}
\keyword{hplot}
