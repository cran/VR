% file MASS/write.matrix.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{write.matrix}
\alias{write.matrix}
\title{
Write a Matrix or Data Frame
}
\description{
Writes a matrix or data frame to a file or the terminal, using column
labels and a layout respecting columns.
}
\usage{
write.matrix(x, file="", sep=" ")
}
\arguments{
\item{x}{
matrix or data frame.
}
\item{file}{
name of output file. The default (\code{""}) is the terminal.
}
\item{sep}{
The separator between columns.
}}
\value{
none
}
\section{Side Effects}{
A formatted file is produced, with column headings (if \code{x} has them)
and columns of data.
}
\seealso{
\code{\link{write.table}}
}
\examples{
\dontrun{
# The function is  defined as
function(x, file = "", sep = " ")
{
   x <- as.matrix(x)
   p <- ncol(x)
   cat(dimnames(x)[[2]], format(t(x)), file = file, 
       sep = c(rep(sep, p - 1), "\n"))
}
}}
\keyword{file}
\keyword{print}
