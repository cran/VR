#-*- R -*-

# Chapter 2   The S Language

library(MASS)
#postscript(file="ch02.ps", width=8, height=6, pointsize=9)


# 2.1  A concise description of S objects

mydata <- c(2.9, 3.4, 3.4, 3.7, 3.7, 2.8, 2.8, 2.5, 2.4, 2.4)
colours <- c("red", "green", "blue", "white", "black")
x1 <- 25:30
x1
mydata[7]
colours[3]

mydata > 3

names(mydata) <- c('a','b','c','d','e','f','g','h','i','j')
mydata
names(mydata)
mydata["e"]

letters[1:5]
mydata[letters[1:5]]
mydata[mydata > 3]

mydata[-c(3:5)]

mode(mydata)
mode(letters)
mode(sin)
length(mydata)
length(letters)
length(sin)

names(mydata) <- NULL     # remove the names
dim(mydata) <- c(2, 5)
mydata
dim(mydata) <- NULL
matrix(mydata, 2, 5)
matrix(mydata, 2, 5, byrow=T)

Empl <- list(employee="Anna", spouse="Fred", children=3,
               child.ages=c(4,7,9))
Empl$employee
Empl$child.ages[2]
names(Empl) <- letters[1:4]
Empl[3:4]
Empl <- c(Empl, service = 8)
unlist(Empl)
unlist(Empl, use.names=F)

c(list(x = 1:3, a = 3:6), list(y = 8:23, b = c(3, 8, 39)))
c(list(x = 1:3, a = 3:6), list(y = 8:23, b = c(3, 8, 39)), recursive=T)


citizen <- factor(c("uk","us","no","au","uk","us","us"))
citizen
print.default(citizen)
unclass(citizen)

citizen <- factor(c("uk","us","no","au","uk","us","us"),
     levels = c("us", "fr", "no", "au", "uk"))
citizen
table(citizen)

income <- ordered(c("Mid","Hi","Lo","Mid","Lo","Hi","Lo"))
income
as.numeric(income)

inc <- ordered(c("Mid","Hi","Lo","Mid","Lo","Hi","Lo"),
    levels = c("Lo", "Mid", "Hi"))
inc

data(geyser)
erupt <- cut(geyser$duration, breaks = 0:6)
erupt <- ordered(erupt, labels=levels(erupt))
erupt

data(painters)
painters
row.names(painters)
painters[1:5, c(2, 4)]


# 2.2  Calling conventions for functions

args(hist)


# 2.3  Arithmetical expressions

x <- c(10.4, 5.6, 3.1, 6.4, 21.7)
y <- c(x, x)
v <- 2 * x + y + 1
xtrunc <- pmax(0, pmin(1,x))
s3 <- seq(-5, 5, by=0.2)
s3
s4 <- seq(length=51, from=-5, by=0.2)
s4
s5 <- rep(x, times=5)
s5
x <- 1:4          # puts c(1,2,3,4)             into x
x
i <- rep(2, 4)    # puts c(2,2,2,2)             into i
i
y <- rep(x, 2)    # puts c(1,2,3,4,1,2,3,4)     into y
y
z <- rep(x, i)    # puts c(1,1,2,2,3,3,4,4)     into z
z
w <- rep(x, x)    # puts c(1,2,2,3,3,3,4,4,4,4) into w
w

colc <- rep(1:3,rep(8,3));  colc
rowc <- rep(rep(1:4,rep(2,4)), 3); rowc
1 + (ceiling(1:24/8) - 1) %% 3 -> colc; colc
1 + (ceiling(1:24/2) - 1) %% 4 -> rowc; rowc


# 2.4, 2.5 omitted


# 2.6  Character vector operations

paste(c("X","Y"), 1:4)
paste(c("X","Y"), 1:4, sep="")
paste(c("X","Y"), 1:4, sep="", collapse=" + ")
data(state)
substring(state.name[44:50], 1, 4)


# 2.7 omitted


# 2.8  Indexing vectors, matrices and arrays

letters[1:3]
letters[1:3][c(1:3,3:1)]

longitude <- state.center[["x"]]
names(longitude) <- state.name
longitude[c("Hawaii", "Alaska")]

a <- 1:4
a[0]
a[0] <- 10
a


mydata
sort(mydata)

x <- rnorm(100001)
sort(x, partial=50001)[50001]

latitude <- state.center[["y"]]
names(latitude) <- state.name
i <- sort.list(longitude)
cbind(latitude = latitude[i], longitude = longitude[i])

data(shoes)
shoes$B
rank(shoes$B)
rank(round(shoes$B))
sort.list(sort.list(round(shoes$B)))


# 2.9  Input/Output facilities

d <- date()
cat("Today's date is:", substring(d,1,10),
             substring(d,25,28), "\n")
cat(1,2,3,4,5,6, fill=8, labels=letters)
data(iris3)
write(iris3[,1,1], "", 15)
cat(format(iris3[,1,1]), fill=60)

#if(version$major >= 5) showConnections(all=T)

# End of ch02
