\name{parcoord}
\alias{parcoord}
\title{
Parallel Coordinates Plot
}
\description{
Parallel coordinates plot
}
\usage{
parcoord(x, col = 1, lty = 1, \dots)
}
\arguments{
\item{x}{
a matrix or data frame who columns represent variables.
}
\item{col}{
A vector of colours, recycled as necessary for each observation.
}
\item{lty}{
A vector of line types, recycled as necessary for each observation.
}
\item{...}{
Further graphics parameters which are passed to \code{matplot}.
}}
\value{
None.
}
\section{Side Effects}{
a parallel coordinates plots is drawn.
}
\references{
Wegman, E. J. (1990) Hyperdimensional data analysis using parallel
coordinates. 
\emph{Journal of the American Statistical Association,}
\bold{85}, 664--675.
}
\examples{
# The function is currently defined as
function(x, col = 1, lty = 1, \dots)
{
    x <- apply(x, 2, function(x) (x - min(x))/(max(x) - min(x)))
    matplot(1:ncol(x), t(x), type = "l", col = col, lty = lty,
        xlab = "", ylab = "", axes = FALSE, \dots)
    axis(1, at = 1:ncol(x), labels = dimnames(x)[[2]])
    for(i in 1:ncol(x))
        lines(c(i, i), c(0, 1), col = "grey70")
    invisible()
}
}
\keyword{hplot}
