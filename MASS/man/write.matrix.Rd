% file MASS/write.matrix.d
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{write.matrix}
\alias{write.matrix}
\title{
Write a Matrix or Data Frame
}
\description{
Writes a matrix or data frame to a file or the console, using column
labels and a layout respecting columns.
}
\usage{
write.matrix(x, file="", sep=" ", blocksize)
}
\arguments{
\item{x}{
matrix or data frame.
}
\item{file}{
name of output file. The default (\code{""}) is the console.
}
\item{sep}{
The separator between columns.
}
\item{blocksize}{
If supplied and positive, the output is written in blocks of
\code{blocksize} rows.
}}
\details{
Supplying \code{blocksize} is more memory-efficient and enables larger
matrices to be writted, but each block of rows might be formatted
slightly differently.
}
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
\keyword{file}
\keyword{print}
