% file MASS/plot.mca.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{plot.mca}
\alias{plot.mca}
\title{
Plot Method for Objects of Class `mca'
}
\description{
Plot a multiple correspondence analysis.
}
\usage{
plot.mca(x, rows=TRUE, col, cex=rep(par("cex"), 2), \dots)
}
\arguments{
\item{x}{
An object of class \code{"mca"}.
}
\item{rows}{
Should the coordinates for the rows be plotted, or just the vertices
for the levels?
}
\item{col, cex}{
The colours and \code{cex} to be used for the row points and level vertices
respectively.
}
\item{...}{
Additional parameters to \code{plot}.
}}
\seealso{
\code{\link{mca}}, \code{\link{predict.mca}}
}
\examples{
data(farms)
plot(mca(farms, abbrev=TRUE))
}
\keyword{hplot}
\keyword{multivariate}
