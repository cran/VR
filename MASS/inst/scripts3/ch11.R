#-*- R -*-

## Script for the Third Edition of `Modern Applied Statistics with S-PLUS'

# Chapter 11   Multivariate Analysis

# library(mva)
library(MASS) # needed for biplot.princomp
postscript(file="ch11.ps", width=8, height=6, pointsize=9)
options(width=65, digits=5)

data(swiss)
swiss.x <- as.matrix(swiss[,-1])

# 11.1  Graphical methods

data(iris3)
ir <- rbind(iris3[,,1], iris3[,,2], iris3[,,3])
ir.species <- factor(c(rep("s",50), rep("c",50), rep("v",50)))
#if(interactive()) brush(ir)

ir.pca <- princomp(log(ir), cor=T)
ir.pca
summary(ir.pca)
plot(ir.pca)
loadings(ir.pca)
ir.pc <- predict(ir.pca)
eqscplot(ir.pc[,1:2], type="n",
     xlab="first principal component",
     ylab = "second principal component")
text(ir.pc[,1:2], labels = as.character(ir.species))

ir.scal <- cmdscale(dist(ir), k = 2, eig = T)
ir.scal$points[, 2] <- -ir.scal$points[, 2]
eqscplot(ir.scal$points, type="n")
text(ir.scal$points, labels = as.character(ir.species), cex = 0.8)

distp <- dist(ir)
dist2 <- dist(ir.scal$points)
sum((distp - dist2)^2)/sum(distp^2)

ir.sam <- sammon(dist(ir[-143,]))
eqscplot(ir.sam$points, type="n")
text(ir.sam$points, labels = as.character(ir.species[-143]), cex = 0.8)

ir.iso <- isoMDS(dist(ir[-143,]))
eqscplot(ir.iso$points, type="n")
text(ir.iso$points, labels = as.character(ir.species[-143]), cex = 0.8)

data(fgl)
fgl.iso <- isoMDS(dist(as.matrix(fgl[-40, -10])))
eqscplot(fgl.iso$points, type="n", xlab="", ylab="")
# either
for(i in seq(along=levels(fgl$type))) {
  set <- fgl$type[-40] == levels(fgl$type)[i]
  points(fgl.iso$points[set,], pch=18, cex=0.6, col=2+i)}
#key(text=list(levels(fgl$type), col=3:8))
# or
#text(fgl.iso$points, labels = c("F", "N", "V", "C", "T", "H")
#     [fgl$type[-40]], cex=0.6)

data(state)
state <- state.x77[,2:7]
state <- state.x77[,2:7]; row.names(state) <- state.abb
biplot(princomp(state, cor=T), pc.biplot=T, cex = 0.7, ex=0.8)


# 11.2  Cluster analysis

h <- hclust(dist(swiss.x), method="single")
plot(h) # or plclust with library(cluster)
cutree(h, 3)
#plclust( clorder(h, cutree(h, 3) ))

h <- hclust(dist(swiss.x), method="average")
initial <- tapply(swiss.x, list(rep(cutree(h, 3),
   ncol(swiss.x)), col(swiss.x)), mean)
dimnames(initial) <- list(NULL, dimnames(swiss.x)[[2]])
km <- kmeans(swiss.x, initial)
swiss.pca <- princomp(swiss.x)
swiss.pca
swiss.px <- predict(swiss.pca)
dimnames(km$centers)[[2]] <- dimnames(swiss.x)[[2]]
swiss.centers <- predict(swiss.pca, km$centers)
eqscplot(swiss.px[, 1:2], type="n",
  xlab="first principal component",
  ylab="second principal component")
text(swiss.px[,1:2], labels = km$cluster)
points(swiss.centers[,1:2], pch=3, cex=3)
#identify(swiss.px[, 1:2], cex=0.5)

if(F){
h <- mclust(swiss.x, method = "S*")$tree
plclust( clorder(h, cutree(h, 3) ))

h <- mclust(swiss.x, method = "S*", noise=T)
hclass <- mclass(h, 3)
hclass$class
mreloc(hclass, swiss.x, method = "S*", noise=T)
}

library(cluster)
swiss.pam <- pam(swiss.px, 3)
summary(swiss.pam)
eqscplot(swiss.px[, 1:2], type="n",
  xlab="first principal component",
  ylab="second principal component")
text(swiss.px[,1:2], labels = swiss.pam$clustering)
points(swiss.pam$medoid[,1:2], pch=3, cex=3)

fanny(swiss.px, 3)

par(mfrow=c(2,1))
pltree(agnes(swiss.x, method="single"))
pltree(diana(swiss.x))
par(mfrow=c(1,1))


# 11.3 Correspondence analysis

data(caith)
corresp(caith)

dimnames(caith)[[2]] <- c("F", "R", "M", "D", "B")
par(mfcol=c(1,3))
plot(corresp(caith, nf=2)); title("symmetric")
plot(corresp(caith, nf=2), type="rows"); title("rows")
plot(corresp(caith, nf=2), type="col"); title("columns")
par(mfrow=c(1,1))


# 11.4  Discriminant analysis

ir.lda <- lda(log(ir), ir.species)
ir.lda
ir.ld <- predict(ir.lda, dimen=2)$x
eqscplot(ir.ld, type="n", xlab="first linear discriminant",
     ylab = "second linear discriminant")
text(ir.ld, labels = as.character(ir.species[-143]), cex = 0.8)

plot(ir.lda, dimen=1)
plot(ir.lda, type="density", dimen=1)

data(fgl)
par(mfrow=c(1,2), pty="s")
fgl.lda <- lda(type ~ ., fgl)
fgl.ld <- predict(fgl.lda, dimen=2)$x
eqscplot(fgl.ld, type="n", xlab="LD2", ylab="LD1")
# either
for(i in seq(along=levels(fgl$type))) {
  set <- fgl$type[-40] == levels(fgl$type)[i]
  points(fgl.ld[set,], pch=18, cex=0.6, col=2+i)}
#key(text=list(levels(fgl$type), col=3:8))
# or
#text(fgl.ld, labels = c("F", "N", "V", "C", "T", "H")[fgl$type[-40]], cex=0.6)

fgl.rlda <- lda(type ~ ., fgl, method="t")
fgl.rld <- predict(fgl.rlda, dimen=2)$x
eqscplot(fgl.rld, type="n", xlab="LD2", ylab="LD1")
# either
for(i in seq(along=levels(fgl$type))) {
  set <- fgl$type[-40] == levels(fgl$type)[i]
  points(fgl.rld[set,], pch=18, cex=0.6, col=2+i)}
#key(text=list(levels(fgl$type), col=3:8))
# or
#text(fgl.rld, labels = c("F", "N", "V", "C", "T", "H")[fgl$type[-40]], cex=0.6)
par(mfrow=c(1,1), pty="m")


# 11.5 Classification theory

#increase len if you have enough memory.
predplot <- function(object, main="", len=50, ...)
{
   plot(Cushings[,1], Cushings[,2], log="xy", type="n",
     xlab="Tetrahydrocortisone", ylab = "Pregnanetriol", main=main)
   for(il in 1:4) {
     set <- Cushings$Type==levels(Cushings$Type)[il]
     text(Cushings[set, 1], Cushings[set, 2],
          labels=as.character(Cushings$Type[set]), col = 2 + il) }
   xp <- seq(0.6, 4.0, length=len)
   yp <- seq(-3.25, 2.45, length=len)
   cushT <- expand.grid(Tetrahydrocortisone=xp,
     Pregnanetriol=yp)
   Z <- predict(object, cushT, ...); zp <- as.numeric(Z$class)
   zp <- Z$post[,3] - pmax(Z$post[,2], Z$post[,1])
   contour(exp(xp), exp(yp), matrix(zp, len),
     add=T, levels=0, labex=0)
   zp <- Z$post[,1] - pmax(Z$post[,2], Z$post[,3])
   contour(exp(xp), exp(yp), matrix(zp, len),
     add=T, levels=0, labex=0)
   invisible()
}

data(Cushings)
cush <- log(as.matrix(Cushings[, -3]))
tp <- factor(Cushings$Type[1:21])
cush.lda <- lda(cush[1:21,], tp); predplot(cush.lda, "LDA")
cush.qda <- qda(cush[1:21,], tp); predplot(cush.qda, "QDA")
predplot(cush.qda, "QDA (predictive)", method = "predictive")
predplot(cush.qda, "QDA (debiased)", method = "debiased")

library(nnet)
Cf <- data.frame(tp = tp,
  Tetrahydrocortisone = log(Cushings[1:21,1]),
  Pregnanetriol = log(Cushings[1:21,2]) )
cush.multinom <- multinom(tp ~ Tetrahydrocortisone
  + Pregnanetriol, Cf, maxit=250)
xp <- seq(0.6, 4.0, length=100); np <- length(xp)
yp <- seq(-3.25, 2.45, length=100)
cushT <- expand.grid(Tetrahydrocortisone=xp,
    Pregnanetriol=yp)
Z <- predict(cush.multinom, cushT, type="probs")
plot(Cushings[,1], Cushings[,2], log="xy", type="n",
     xlab="Tetrahydrocortisone", ylab = "Pregnanetriol")
for(il in 1:4) {
  set <- Cushings$Type==levels(Cushings$Type)[il]
  text(Cushings[set, 1], Cushings[set, 2],
       labels=as.character(Cushings$Type[set]), col = 2 + il) }
zp <- Z[,3] - pmax(Z[,2], Z[,1])
contour(exp(xp), exp(yp), matrix(zp, np),
   add=T, levels=0, labex=0)
zp <- Z[,1] - pmax(Z[,2], Z[,3])
contour(exp(xp), exp(yp), matrix(zp, np),
   add=T, levels=0, labex=0)

library(tree)
cush.tr <- tree(tp ~ Tetrahydrocortisone + Pregnanetriol, Cf)
plot(cush[,1], cush[,2], type="n",
     xlab="Tetrahydrocortisone", ylab = "Pregnanetriol")
for(il in 1:4) {
  set <- Cushings$Type==levels(Cushings$Type)[il]
  text(cush[set, 1], cush[set, 2],
       labels= as.character(Cushings$Type[set]), col = 2 + il) }
par(cex=1.5); partition.tree(cush.tr, add=T); par(cex=1)


# 11.6  Other classification methods

# neural networks

library(nnet)
pltnn <- function(main, ...) {
   plot(Cushings[,1], Cushings[,2], log="xy", type="n",
   xlab="Tetrahydrocortisone", ylab = "Pregnanetriol", main=main, ...)
   for(il in 1:4) {
       set <- Cushings$Type==levels(Cushings$Type)[il]
       text(Cushings[set, 1], Cushings[set, 2],
          as.character(Cushings$Type[set]), col = 2 + il) }
}

plt.bndry <- function(size=0, decay=0, ...)
{
   cush.nn <- nnet(cush, tpi, skip=T, softmax=T, size=size,
      decay=decay, maxit=1000)
   invisible(b1(predict(cush.nn, cushT), ...))
}

b1 <- function(Z, ...)
{
   zp <- Z[,3] - pmax(Z[,2], Z[,1])
   contour(exp(xp), exp(yp), matrix(zp, np),
      add=T, levels=0, labex=0, ...)
   zp <- Z[,1] - pmax(Z[,3], Z[,2])
   contour(exp(xp), exp(yp), matrix(zp, np),
      add=T, levels=0, labex=0, ...)
}

cush <- cush[1:21,]; tpi <- class.ind(tp)
# functions pltnn and plt.bndry given in the scripts
par(mfrow=c(2,2))
pltnn("Size = 2")
set.seed(1); plt.bndry(size=2, col=2)
set.seed(3); plt.bndry(size=2, col=3); plt.bndry(size=2, col=4)

pltnn("Size = 2, lambda = 0.001")
set.seed(1); plt.bndry(size=2, decay=0.001, col=2)
set.seed(2); plt.bndry(size=0, decay=0.001, col=4)

pltnn("Size = 2, lambda = 0.01")
set.seed(1); plt.bndry(size=2, decay=0.01, col=2)
set.seed(2); plt.bndry(size=2, decay=0.01, col=4)

pltnn("Size = 5, 20  lambda = 0.01")
set.seed(2); plt.bndry(size=5, decay=0.01, col=1)
set.seed(2); plt.bndry(size=20, decay=0.01, col=2)

# functions pltnn and b1 are in the scripts
pltnn("Many local maxima")
Z <- matrix(0, nrow(cushT), ncol(tpi))
for(iter in 1:20) {
   set.seed(iter)
   cush.nn <- nnet(cush, tpi,  skip=T, softmax=T, size=3,
       decay=0.01, maxit=1000, trace=F)
   Z <- Z + predict(cush.nn, cushT)
# In 5.x replace $ by @ in next line.
   cat("final value", format(round(cush.nn$value,3)), "\n")
   b1(predict(cush.nn, cushT), col=2, lwd=0.5)
}
pltnn("Averaged")
b1(Z, lwd=3)

# Non-parametric rules

library(class)
par(pty="s", mfrow=c(1,2))
plot(Cushings[,1], Cushings[,2], log="xy", type="n",
     xlab = "Tetrahydrocortisone", ylab = "Pregnanetriol",
     main = "1-NN")
for(il in 1:4) {
  set <- Cushings$Type==levels(Cushings$Type)[il]
  text(Cushings[set, 1], Cushings[set, 2],
       as.character(Cushings$Type[set]), col = 2 + il) }
Z <- knn(scale(cush, F, c(3.4, 5.7)),
         scale(cushT, F, c(3.4, 5.7)), tp)
contour(exp(xp), exp(yp), matrix(as.numeric(Z=="a"), np),
      add=T, levels=0.5, labex=0)
contour(exp(xp), exp(yp), matrix(as.numeric(Z=="c"), np),
      add=T, levels=0.5, labex=0)
plot(Cushings[,1], Cushings[,2], log="xy", type="n",
     xlab="Tetrahydrocortisone", ylab = "Pregnanetriol",
     main = "3-NN")
for(il in 1:4) {
  set <- Cushings$Type==levels(Cushings$Type)[il]
  text(Cushings[set, 1], Cushings[set, 2],
       as.character(Cushings$Type[set]), col = 2 + il) }
Z <- knn(scale(cush, F, c(3.4, 5.7)),
         scale(cushT, F, c(3.4, 5.7)), tp, k=3)
contour(exp(xp), exp(yp), matrix(as.numeric(Z=="a"), np),
      add=T, levels=0.5, labex=0)
contour(exp(xp), exp(yp), matrix(as.numeric(Z=="c"), np),
      add=T, levels=0.5, labex=0)
par(pty="m", mfrow=c(1,1))


# 11.7 Two extended examples

# Leptograpsus variegatis crabs

data(crabs)
lcrabs <- log(crabs[,4:8])
crabs.grp <- factor(c("B", "b", "O", "o")[rep(1:4, rep(50,4))])
lcrabs.pca <- princomp(lcrabs)
lcrabs.pc <- predict(lcrabs.pca)
dimnames(lcrabs.pc) <- list(NULL, paste("PC", 1:5, sep=""))
lcrabs.pca
loadings(lcrabs.pca)

if(F) {
library(xgobi)
xgobi(lcrabs, colors=c("SkyBlue", "SlateBlue", "Orange",
      "Red")[rep(1:4, rep(50, 4))])
xgobi(lcrabs, glyphs=12 + 5*rep(0:3, rep(50, 4)))
}

cr.scale <- 0.5 * log(crabs$CL * crabs$CW)
slcrabs <- lcrabs - cr.scale
cr.means <- matrix(0, 2, 5)
cr.means[1,] <- apply(slcrabs[crabs$sex=="F",], 2, mean)
cr.means[2,] <- apply(slcrabs[crabs$sex=="M",], 2, mean)
dslcrabs <- slcrabs - cr.means[as.numeric(crabs$sex),]
lcrabs.sam <- sammon(dist(dslcrabs))
eqscplot(lcrabs.sam$points, type="n", xlab="", ylab="")
text(lcrabs.sam$points, labels = as.character(crabs.grp))

if(F) {
crabs.h <- cutree(hclust(dist(dslcrabs)),2)
table(crabs$sp, crabs.h)
cr.means[1,] <- apply(dslcrabs[crabs.h==1,], 2, mean)
cr.means[2,] <- apply(dslcrabs[crabs.h==2,], 2, mean)
crabs.km <- kmeans(dslcrabs, cr.means)
table(crabs$sp, crabs.km$cluster)
}

dcrabs.lda <- lda(crabs$sex ~ FL + RW + CL + CW, lcrabs)
dcrabs.lda
dcrabs.pred <- predict(dcrabs.lda)
table(crabs$sex, dcrabs.pred$class)

dcrabs.lda4 <- lda(crabs.grp ~ FL + RW + CL + CW, lcrabs)
dcrabs.lda4
dcrabs.pr4 <- predict(dcrabs.lda4, dimen=2)
dcrabs.pr2 <- dcrabs.pr4$post[, c("B","O")] %*% c(1,1)
table(crabs$sex, dcrabs.pr2 > 0.5)

cr.t <- dcrabs.pr4$x[,1:2]
eqscplot(cr.t, type="n", xlab="First LD", ylab="Second LD")
text(cr.t, labels = as.character(crabs.grp))
perp <- function(x, y) {
   m <- (x+y)/2
   s <- - (x[1] - y[1])/(x[2] - y[2])
   abline(c(m[2] - s*m[1], s))
   invisible()
}
# For 5.x replace $means by @means
cr.m <- lda(cr.t, crabs$sex)$means
points(cr.m, pch=3, mkh=0.3)
perp(cr.m[1,], cr.m[2,])

cr.lda <- lda(cr.t, crabs.grp)
x <- seq(-6, 6, 0.25)
y <- seq(-2, 2, 0.25)
Xcon <- matrix(c(rep(x,length(y)),
               rep(y, rep(length(x),length(y)))),,2)
cr.pr <- predict(cr.lda, Xcon)$post[, c("B","O")] %*% c(1,1)
contour(x, y, matrix(cr.pr, length(x), length(y)),
  levels=0.5, labex=0, add=T, lty=3)

for(i in c("O", "o",  "B", "b"))
    print(var(lcrabs[crabs.grp==i, ]))

# Forensic glass

data(fgl)
set.seed(123); rand <- sample (10, 214, replace=T)
con <- function(x,y)
{
  tab <- table(x,y)
  print(tab)
  diag(tab) <- 0
  cat("error rate = ", round(100*sum(tab)/length(x),2),"%\n")
  invisible()
}

CVtest <- function(fitfn, predfn, ...)
{
  res <- fgl$type
  for (i in sort(unique(rand))) {
    cat("fold ",i,"\n", sep="")
    learn <- fitfn(rand != i, ...)
    res[rand == i] <- predfn(learn, rand==i)
  }
  res
}

res.multinom <- CVtest(
  function(x, ...) multinom(type ~ ., fgl[x,], ...),
  function(obj, x) predict(obj, fgl[x, ],type="class"),
  maxit=1000, trace=F )
con(fgl$type, res.multinom)

res.lda <- CVtest(
  function(x, ...) lda(type ~ ., fgl[x, ], ...),
  function(obj, x) predict(obj, fgl[x, ])$class )
con(fgl$type, res.lda)

library(class)
fgl0 <- fgl[ ,-10] # drop type
{ res <- fgl$type
 for (i in sort(unique(rand))) {
    cat("fold ",i,"\n", sep="")
    sub <- rand == i
    res[sub] <- knn(fgl0[!sub, ], fgl0[sub,], fgl$type[!sub], k=1)
 }
 res } -> res.knn1
con(fgl$type, res.knn1)

res.lb <- knn(fgl0, fgl0, fgl$type, k=3, prob=T, use.all=F)
table(attr(res.lb, "prob"))

library(rpart) ## rpart 3.0 or later
res.rpart <- CVtest(
  function(x, ...) {
    tr <- rpart(type ~ ., fgl[x,], ...)
    cp <- tr$cptable
    r <- cp[, 4] + cp[, 5]
    rmin <- min(seq(along=r)[cp[, 4] < min(r)])
    cp0 <- cp[rmin, 1]
    cat("size chosen was", cp[rmin, 2] + 1, "\n")
    prune(tr, cp=1.01*cp0)
  },
  function(obj, x)
     predict(obj, fgl[x, ], type="class"),
#    levels(fgl$type)[apply(predict(obj, fgl[x, ]), 1, which.is.max)],
  cp = 0.001
)
con(fgl$type, res.rpart)

fgl1 <- lapply(fgl[, 1:9], function(x) {r <- range(x); (x-r[1])/diff(r)})
fgl1 <- data.frame(fgl1, type=fgl$type)

# stepAIC does not currently work with new-style classes like multinom
res.mult3 <- CVtest(
  function(xsamp, ...) {
    assign("xsamp", xsamp, envir=.GlobalEnv)
    obj <- multinom(type ~ ., fgl1[xsamp,], trace=F, ...)
    stepAIC(obj)
  },
  function(obj, x) predict(obj, fgl1[x, ],type="class"),
  maxit=1000, decay=1e-3)
con(fgl$type, res.mult3)

CVnn2 <- function(formula, data,
                  size = rep(6,2), lambda = c(0.001, 0.01),
                  nreps = 1, nifold = 5, verbose = 99, ...)
{
  CVnn1 <- function(formula, data, nreps=1, ri, verbose,  ...)
  {
    truth <- data[,deparse(formula[[2]])]
    res <-  matrix(0, nrow(data), length(levels(truth)))
    if(verbose > 20) cat("  inner fold")
    for (i in sort(unique(ri))) {
      if(verbose > 20) cat(" ", i,  sep="")
      for(rep in 1:nreps) {
        learn <- nnet(formula, data[ri !=i,], trace=F, ...)
        res[ri == i,] <- res[ri == i,] +
          predict(learn, data[ri == i,])
      }
   }
    if(verbose > 20) cat("\n")
    sum(as.numeric(truth) != max.col(res/nreps))
  }
  truth <- data[,deparse(formula[[2]])]
  res <-  matrix(0, nrow(data), length(levels(truth)))
  choice <- numeric(length(lambda))
  for (i in sort(unique(rand))) {
    if(verbose > 0) cat("fold ", i,"\n", sep="")
    ri <- sample(nifold, sum(rand!=i), replace=T)
    for(j in seq(along=lambda)) {
      if(verbose > 10)
        cat("  size =", size[j], "decay =", lambda[j], "\n")
      choice[j] <- CVnn1(formula, data[rand != i,], nreps=nreps,
                          ri=ri, size=size[j], decay=lambda[j],
                          verbose=verbose, ...)
    }
    decay <- lambda[which.is.max(-choice)]
    csize <- size[which.is.max(-choice)]
    if(verbose > 5) cat("  #errors:", choice, "  ")
    if(verbose > 1) cat("chosen size = ", csize,
                        " decay = ", decay, "\n", sep="")
    for(rep in 1:nreps) {
      learn <- nnet(formula, data[rand != i,], trace=F,
                    size=csize, decay=decay, ...)
      res[rand == i,] <- res[rand == i,] +
          predict(learn, data[rand == i,])
    }
  }
  factor(levels(truth)[max.col(res/nreps)], levels = levels(truth))
}

if(F) { # only run this if you have time to wait
res.nn2 <- CVnn2(type ~ ., fgl1, skip=T, maxit=500, nreps=10)
con(fgl$type, res.nn2)
}

cd0 <- lvqinit(fgl0, fgl$type, prior=rep(1,6)/6,k=3)
cd1 <- olvq1(fgl0, fgl$type, cd0)
con(fgl$type, lvqtest(cd1, fgl0))

CV.lvq <- function()
{
 res <- fgl$type
 for(i in sort(unique(rand))) {
   cat("doing fold",i,"\n")
   cd0 <- lvqinit(fgl0[rand != i,], fgl$type[rand != i],
                  prior=rep(1,6)/6, k=3)
   cd1 <- olvq1(fgl0[rand != i,], fgl$type[rand != i], cd0)
   cd1 <- lvq3(fgl0[rand != i,], fgl$type[rand != i],
               cd1, niter=10000)
   res[rand == i] <- lvqtest(cd1, fgl0[rand == i,])
 }
 res
}
con(fgl$type, CV.lvq())

# Try Mahalanobis distance
# library(mva)
fgl0 <- scale(princomp(fgl[,-10])$scores)
con(fgl$type, CV.lvq())


# 11.8 Calibration plots

CVprobs <- function(fitfn, predfn, ...)
{
 res <- matrix(, 214, 6)
 for (i in sort(unique(rand))) {
    cat("fold ",i,"\n", sep="")
    learn <- fitfn(rand != i, ...)
    res[rand == i,] <- predfn(learn, rand==i)
 }
 res
}
probs.multinom <- CVprobs(
  function(x, ...) multinom(type ~ ., fgl[x,], ...),
  function(obj, x) predict(obj, fgl[x, ],type="probs"),
  maxit=1000, trace=F )

# library(modreg)
probs.yes <- as.vector(class.ind(fgl$type))
probs <- as.vector(probs.multinom)
par(pty="s")
plot(c(0,1), c(0,1), type="n", xlab="predicted probability",
  ylab="", xaxs="i", yaxs="i", las=1)
rug(probs[probs.yes==0], 0.02, side=1, lwd=0.5)
rug(probs[probs.yes==1], 0.02, side=3, lwd=0.5)
abline(0,1)
newp <- seq(0, 1, length=100)
lines(newp, predict(loess(probs.yes ~ probs, span=1), newp))

# End of ch11
