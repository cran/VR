% file nnet/nnet.Hess.d copyright (C) 1994-2000 W. N. Venables and B. D. Ripley
%
\name{nnet.Hess}
\alias{nnet.Hess}
\title{
Evaluates Hessian for a Neural Network
}
\description{
Evaluates the Hessian (matrix of second derivatives) of the specified
neural network. Normally called via argument \code{Hess=TRUE} to \code{nnet} or via
\code{vcov.multinom}.
}
\usage{
nnet.Hess(net, x, y, weights)
}
\arguments{
\item{net}{
object of class \code{nnet} as returned by \code{nnet}.
}
\item{x}{
training data.
}
\item{y}{
classes for training data.
}
\item{weights}{
the (case) weights used in the \code{nnet} fit.
}}
\value{
square symmetric matrix of the Hessian evaluated at the weights stored 
in the net.
}
\seealso{
\code{\link{nnet}}, \code{\link{predict.nnet}}
}
\examples{
data(iris3)
# use half the iris data
ir <- rbind(iris3[,,1],iris3[,,2],iris3[,,3])
targets <- matrix(c(rep(c(1,0,0),50), rep(c(0,1,0),50), rep(c(0,0,1),50)),
150, 3, byrow=TRUE)
samp <- c(sample(1:50,25), sample(51:100,25), sample(101:150,25))
ir1 <- nnet(ir[samp,], targets[samp,],size=2, rang=0.1, decay=5e-4, maxit=200)
eigen(nnet.Hess(ir1, ir[samp,], targets[samp,]), T)$values
}
\keyword{neural}
