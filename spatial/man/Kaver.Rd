% file spatial/Kaver.d copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{Kaver}
\alias{Kaver}
\title{
Average K-functions from Simulations
}
\description{
Forms the average of a series of (usually simulated) K-functions.
}
\usage{
Kaver(fs, nsim, \dots)
}
\arguments{
\item{fs}{
full scale for K-fn
}
\item{nsim}{
number of simulations
}
\item{...}{
arguments to simulate one point process object
}}
\value{
list with components \code{x} and \code{y} of the average K-fn on L-scale.
}
\seealso{
\code{\link{Kfn}}, \code{\link{Kenvl}}
}
\examples{
towns <- ppinit("towns.dat")
par(pty="s")
plot(Kfn(towns, 40), type="b")
plot(Kfn(towns, 10), type="b", xlab="distance", ylab="L(t)")
for(i in 1:10) lines(Kfn(Psim(69), 10))
lims <- Kenvl(10,100,Psim(69))
lines(lims$x,lims$l, lty=2, col="green")
lines(lims$x,lims$u, lty=2, col="green")
lines(Kaver(10,25,Strauss(69,0.5,3.5)),  col="red")
}
\keyword{spatial}
