% file MASS/kde2d.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{kde2d}
\alias{kde2d}
\title{
Two-Dimensional Kernel Density Estimation
}
\description{
Two-dimensional kernel density estimation with an axis-aligned
bivariate normal kernel, evaluated on a square grid.
}
\usage{
kde2d(x, y, h, n=25, lims=c(range(x), range(y)))
}
\arguments{
\item{x}{
x coordinate of data
}
\item{y}{
y coordinate of data
}
\item{h}{
vector of bandwidths for x and y directions. Defaults to
normal reference bandwidth.
}
\item{n}{
Number of grid points in each direction.
}
\item{lims}{
The limits of the rectangle covered by the grid as \code{c(xl, xu, yl, yu)}.
}}
\value{
A list of three components.

\item{x,}{
The x and y coordinates of the grid points, vectors of length \code{n}.
}
\item{z}{
An \code{n}  x  \code{n} matrix of the evaluated density.
}}
keyword{dplot}
% Converted by Sd2Rd version 0.3-2.
