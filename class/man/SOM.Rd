\name{SOM}
\alias{SOM}
\title{
Self-Organizing Maps
}
\description{
Kohonen's Self-Organizing Maps are a crude form of multidimensional scaling.
}
\usage{
SOM(data, grid = somgrid(), rlen = 10000, alpha, radii, init)
}
\arguments{
\item{data}{
a matrix or data frame of observations, scaled so that Euclidean
distance is appropriate.
}
\item{grid}{
A grid for the representatives.
}
\item{rlen}{
the number of updates: used only in the defaults for \code{alpha} and \code{radii}. 
}
\item{alpha}{
the amount of change: one update is done for each element of \code{alpha}.
Default is to decline linearly from 0.05 to 0 over \code{rlen} updates.
}
\item{radii}{
the radii of the neighbourhood to be used for each update: must be the
same length as \code{alpha}.  Default is to decline linearly from 4 to 1
over \code{rlen} updates.
}
\item{init}{
the initial representatives.  If missing, chosen (without replacement)
randomly from \code{data}.
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
\code{alpha} and \code{radii} can also be lists, in which case each component is
used in turn, allowing two- or more phase training.
}
\seealso{
\code{\link{batchSOM}}
}
\references{
Kohonen, T. (1995)
\emph{Self-Organizing Maps.}
Springer-Verlag

Kohonen, T., Hynninen, J., Kangas, J. and Laaksonen, J. (1996)
SOM PAK: The self-organizing map program package.
Laboratory of Computer and Information Science, Helsinki University
of Technology, Technical Report A31
}
\examples{
library(MASS)
data(crabs)
library(mva) # for dist

lcrabs <- log(crabs[, 4:8])
crabs.grp <- factor(c("B", "b", "O", "o")[rep(1:4, rep(50,4))])
gr <- somgrid(topo = "hexagonal")
crabs.som <- SOM(lcrabs, gr)
plot(crabs.som)

## 2-phase training
crabs.som2 <- SOM(lcrabs, gr,
    alpha = list(seq(0.05, 0, len = 1e4), seq(0.02, 0, len = 1e5)),
    radii = list(seq(8, 1, len = 1e4), seq(4, 1, len = 1e5)))
plot(crabs.som2)
}
\keyword{classif}
