\name{batchSOM}
\alias{batchSOM}
\alias{somgrid}
\alias{plot.somgrid}
\alias{plot.SOM}
\title{
Self-Organizing Maps
}
\description{
Kohonen's Self-Organizing Maps are a crude form of multidimensional scaling.
}
\usage{
batchSOM(data, grid = somgrid(), radii, init)
somgrid(xdim = 8, ydim = 6, topo = c("rectangular", "hexagonal"))
}
\arguments{
\item{data}{
a matrix or data frame of observations, scaled so that Euclidean
distance is appropriate.
}
\item{grid}{
A grid for the representatives.
}
\item{radii}{
the radii of the neighbourhood to be used for each pass: one pass is
run for each element of \code{radii}.
}
\item{init}{
the initial representatives.  If missing, chosen (without replacement)
randomly from \code{data}.
}
\item{xdim, ydim}{
dimensions of the grid
}
\item{topo}{
the topology of the grid.
}}
\value{
an object of class \code{"SOM"} with components

\item{grid}{
the grid, an object of class \code{"somgrid"}.
}
\item{codes}{
a matrix of representatives.
}}
\details{
The batch SOM algorithm of Kohonen(1995, section 3.14) is used.
}
\references{
Kohonen, T. (1995)
\emph{Self-Organizing Maps.}
Springer-Verlag
}
\examples{
library(MASS)
data(crabs)
library(mva) # for dist
library(class)
lcrabs <- log(crabs[, 4:8])
crabs.grp <- factor(c("B", "b", "O", "o")[rep(1:4, rep(50,4))])
gr <- somgrid(topo = "hexagonal")
crabs.som <- batchSOM(lcrabs, gr, c(4, 4, 2, 2, 1, 1, 1, 0, 0))
plot(crabs.som)

bins <- as.numeric(knn1(crabs.som$code, lcrabs, 0:47))
plot(crabs.som$grid, type = "n")
symbols(crabs.som$grid$pts[, 1], crabs.som$grid$pts[, 2], 
        circles = rep(0.4, 48), inches = FALSE, add = TRUE)
text(crabs.som$grid$pts[bins, ] + rnorm(400, 0, 0.1), 
     as.character(crabs.grp))
}
\keyword{classif}
