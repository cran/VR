#-*- R -*-

# Chapter 4   Programming in S

library(MASS)


# 4.1  Control structures

yp <- rpois(50, lam=1)   # full Poisson sample of size 50
table(yp)
y <- yp[yp > 0]          # truncate the zeros
ybar <- mean(y); ybar
lam <- ybar
it <- 0                  # iteration count
del <- 1                 # iterative adjustment
while (abs(del) > 0.0001 && (it <- it + 1) < 10) {
    del <- (lam - ybar*(1 - exp(-lam)))/(1 - ybar*exp(-lam))
    lam <- lam - del
    cat(it, lam, "\n")
    }


# 4.2  More on character strings

data(state)
as.vector(abbreviate(state.name[44:50]))
as.vector(abbreviate(state.name[44:50], use.classes=F))

if(F) {
# S-PLUS 3.4 and later have regexpr
regexpr("na$", state.name)
state.name[regexpr("na$", state.name)> 0]

if(version$major >= 5) {
  print(regMatch(state.name, "na$"))
  print(regMatchPos(state.name, "na$"))
  print(state.name[regularExpression("na$")])
  print(substring(state.name, "^[A-Za-z]+"))
}
}

# 4.3  Matrix operations

p <- dbinom(0:4, size=4, prob=1/3)  # an example prob vector
CC <- -(p %o% p)
diag(CC) <- p + diag(CC)
structure(3^8 * CC, dimnames=list(0:4, 0:4))  # convenience


# 4.4  Vectorized calculations and loop avoidance functions

data(iris3)
apply(iris3, c(2,3), mean)
apply(iris3, c(2,3), mean, trim=0.1)
apply(iris3, 2, mean)
ir.var <- apply(iris3, 3, var)
ir.var <- array(ir.var, dim = dim(iris3)[c(2,2,3)],
                dimnames = dimnames(iris3)[c(2,2,3)])

matrix(rep(1/50, 50) %*% matrix(iris3, nrow = 50), nrow = 4,
         dimnames = dimnames(iris3)[-1])

ir.means <- apply(iris3, c(2,3), mean)
sweep(iris3, c(2,3), ir.means)
log(sweep(iris3, c(2,3), ir.means, "/"))

data(quine)
attach(quine)
table(Age)
table(Sex, Age)

#tab <- crosstabs(~Sex + Age, quine)
#print.default(tab)
## R equivalent (from 1.2.0)
xtabs(~Sex + Age, quine)

tapply(Days, Age, mean)
tapply(Days, Age, mean, trim = 0.1)
tapply(Days, list(Sex, Age), mean)
tapply(Days, list(Sex, Age), function(x) sqrt(var(x)/length(x)))


conv1 <- function(a, b) {
   ab <- outer(a, b)
   unlist(lapply(split(ab, row(ab) + col(ab)), sum))
}
conv2 <- function(a, b) {
    ab <- outer(a, b)
    tapply(ab, row(ab) + col(ab), sum)
}

a <- runif(1000)
b <- runif(100)
print(system.time(conv1(a, b)))
print(system.time(conv2(a, b)))

Letters <- c(LETTERS, letters)
Letters[sapply(Letters, function(xx) exists(xx))]
Letters[sapply(Letters, function(x) exists(x))]
Letters[sapply(Letters, exists)]

sapply(list(quine, quine), function(x) dim(x))

sapply(quine, is.factor)
quineFO <- quine[,sapply(quine, is.factor)]
tab <- do.call("table", quineFO) # or just table(quineFO) in R
tab

QuineF <- expand.grid(lapply(quineFO, levels))
QuineF <- do.call("expand.grid", lapply(quineFO, levels))
QuineF$Freq <- as.vector(tab)
QuineF

data(crabs)
aggregate(crabs[, 4:8], list(sp=crabs$sp, sex=crabs$sex), median)
## by and merge need R >= 1.0.0
by(crabs[,4:8], list(species=crabs$sp, sex=crabs$sex), summary)
data(Animals, mammals)
merge(Animals, mammals, by="row.names")


# End of ch04
