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
eqscplot(x, y, tol=0.04, \dots)
}
\arguments{
\item{x}{
vector of x values, or a 2-column matrix, or a list with components
\code{x} and \code{y}
}
\item{y}{
vector of y values
}
\item{tol}{
proportion of "white space" at the margins of plot
}
\item{...}{
further arguments for \code{plot}
}}
\value{
none
}
\section{Side Effects}{
performs the plot.
}
\details{
Limits for the x and y axes are chosen so that they include the data. One of
the sets of limits is then stretched from the midpoint to make the
units equal. Finally both are stretched by \code{1 + tol} to move points away from
the axes, and the points plotted. 
}
\seealso{
\code{\link{plot}}
}
\keyword{hplot}
