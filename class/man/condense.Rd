% file class/condense.d copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{condense}
\alias{condense}
\title{
Condense training set for k-NN classifier
}
\description{
Condense training set for k-NN classifier
}
\usage{
condense(train, class, store)
}
\arguments{
\item{train}{
matrix for training set
}
\item{class}{
vector of classifications for test set
}
\item{store}{
initial store set. Default one randomly chosen element of the set.
}}
\details{
The store set is used to 1-NN classify the rest, and misclassified
patterns are added to the store set. The whole set is checked until
no additions occur.
}
\value{
index vector of cases to be retained (the final store set).
}
\references{
P. A. Devijver and J. Kittler (1982)
\emph{Pattern Recognition. A Statistical Approach.}
Prentice-Hall, pp. 119--121.
}
\seealso{
\code{\link{reduce.nn}}, \code{\link{multiedit}}
}
\examples{
data(iris3)
train <- rbind(iris3[1:25,,1],iris3[1:25,,2],iris3[1:25,,3])
test <- rbind(iris3[26:50,,1],iris3[26:50,,2],iris3[26:50,,3])
cl <- factor(c(rep("s",25),rep("c",25), rep("v",25)))
keep <- condense(train, cl)
knn(train[keep,], test, cl[keep])
keep2 <- reduce.nn(train, keep, cl)
knn(train[keep2,], test, cl[keep2])
}
\keyword{classif}
