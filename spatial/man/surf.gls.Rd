% file spatial/surf.gls.d copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{surf.gls}
\alias{surf.gls}
\title{
Fits a Trend Surface by Generalized Least-squares
}
\description{
Fits a trend surface by generalized least-squares.
}
\usage{
surf.gls(np, covmod, x, y, z, nx = 1000, \dots)
}
\arguments{
\item{np}{
degree of polynomial surface
}
\item{covmod}{
function to evaluate covariance or correlation function
}
\item{x}{
x coordinates or a data frame with columns \code{x}, \code{y}, \code{z}
}
\item{y}{
y coordinates
}
\item{z}{
z coordinates. Will supersede \code{x$z}
}
\item{nx}{
Number of bins for table of the covariance. Increasing adds accuracy,
and increases size of the object.
}
\item{...}{
parameters for \code{covmod}
}}
\value{
list with components

\item{beta}{
the coefficients
}
\item{x}{
}
\item{y}{
}
\item{z}{


and others for internal use only.
}}
\seealso{
\code{\link{trmat}}, \code{\link{surf.ls}}, \code{\link{prmat}}, \code{\link{semat}}, \code{\link{expcov}}, \code{\link{gaucov}}, \code{\link{sphercov}}
}
\examples{
data(topo)
topo.kr <- surf.gls(2, expcov, topo, d=0.7)
trsurf <- trmat(topo.kr, 0, 6.5, 0, 6.5, 50)
contour(trsurf)
prsurf <- prmat(topo.kr, 0, 6.5, 0, 6.5, 50)
contour(prsurf, levels=seq(700, 925, 25))
sesurf <- semat(topo.kr, 0, 6.5, 0, 6.5, 30)
contour(sesurf, levels=c(22,25))
}
\keyword{spatial}
