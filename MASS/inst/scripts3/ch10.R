#-*- R -*-

## Script for the Third Edition of `Modern Applied Statistics with S-PLUS'

# Chapter 10   Tree-based Methods

library(MASS)
library(tree)
postscript(file="ch10.ps", width=8, height=6, pointsize=9)

#set.seed(123)
cpus.samp <- sample(1:209, 100)
data(cpus)

# 10.2  Implementation in rpart

library(rpart)
#set.seed(123)
cpus.rp <- rpart(log10(perf) ~ ., cpus[ ,2:8], cp=1e-3)
cpus.rp  # gives a large tree not show here.
print(cpus.rp, cp=0.01)

plot(cpus.rp, uniform=T); text(cpus.rp, digits=3)

post(cpus.rp, filename="CpusTree.eps", horizontal=F, pointsize=8)

data(iris3)
ird <- data.frame(rbind(iris3[,,1], iris3[,,2],iris3[,,3]),
           Species=c(rep("s",50), rep("c",50), rep("v",50)))
ir.rp <- rpart(Species ~ ., data=ird, method="class", cp=1e-3)
ir.rp

printcp(cpus.rp)

print(cpus.rp, cp=0.006, digits=3)
cpus.rp1 <- prune(cpus.rp, cp=0.006)
plot(cpus.rp1, branch=0.4, uniform=T)
text(cpus.rp1, digits=3)

printcp(ir.rp)

summary(ir.rp)

ir.rp1 <- rpart(Species ~ ., ird, cp=0, minsplit=5, maxsurrogate=0)
printcp(ir.rp1)
print(ir.rp1, cp=0.015)

plot(cpus.rp, branch=0.6, compress=T, uniform=T)
text(cpus.rp, digits=3, all=T, use.n=T)

data(fgl)
#set.seed(123)
fgl.rp <- rpart(type ~ ., fgl, cp=0.001)
plotcp(fgl.rp)
printcp(fgl.rp)
print(fgl.rp, cp=0.02)
fgl.rp2 <- prune(fgl.rp, cp=0.02)
plot(fgl.rp2, uniform=T); text(fgl.rp2, use.n=T)


# 10.3 Implementation in tree

cpus.tr <- tree(perf ~ syct+mmin+mmax+cach+chmin+chmax, cpus)
summary(cpus.tr)
print(cpus.tr)

par(mfrow=c(1,2))
plot(cpus.tr, type="u");  text(cpus.tr, srt=90)

cpus.ltr <- update(cpus.tr, log10(perf) ~ .)
summary(cpus.ltr)
plot(cpus.ltr, type="u");  text(cpus.ltr, srt=90)
par(mfrow = c(1,1))

ir.tr <- tree(Species ~., ird)
summary(ir.tr)
ir.tr
plot(ir.tr)
text(ir.tr, all=T)
ir.tr1 <- snip.tree(ir.tr, nodes = c(12, 7))
ir.tr1
summary(ir.tr1)

par(pty="s")
plot(ird[, 3],ird[, 4], type="n",
  xlab="petal length", ylab="petal width")
text(ird[, 3], ird[, 4], labels=as.character(ird$Species))
par(cex=2)
partition.tree(ir.tr1, add=T)
par(cex=1)

par(pty="m")

fgl.tr <- tree(type ~ ., fgl)
summary(fgl.tr)
plot(fgl.tr);  text(fgl.tr, all=T, cex=0.5)
fgl.tr1 <- snip.tree(fgl.tr, node=c(11, 53, 105, 108, 31))
tree.screens()
plot(fgl.tr1)
tile.tree(fgl.tr1, fgl$type)
close.screen(all = T)

data(shuttle)
shuttle.tr <- tree(use ~ ., shuttle, subset=1:253,
                   mindev=1e-6, minsize=2)
shuttle.tr
#post.tree(shuttle.tr)
shuttle1 <- shuttle[254:256, ]  # 3 missing cases
predict(shuttle.tr, shuttle1)



par(mfrow=c(1,2), pty="s")
plot(prune.tree(cpus.ltr))
cpus.ltr1 <- prune.tree(cpus.ltr, best=8)
plot(cpus.ltr1);   text(cpus.ltr1)
sqrt(sum((log10(cpus[-cpus.samp, "perf"]) -
          predict(cpus.ltr1, cpus[-cpus.samp,]))^2)/109)

summary(prune.tree(fgl.tr, k=10))

#set.seed(123)
plot(cv.tree(cpus.ltr,, prune.tree))
#post.tree(prune.tree(cpus.ltr, best=4))

#set.seed(123)
fgl.cv <- cv.tree(fgl.tr,, prune.tree)
for(i in 2:5)  fgl.cv$dev <- fgl.cv$dev +
   cv.tree(fgl.tr,, prune.tree)$dev
fgl.cv$dev <- fgl.cv$dev/5
plot(fgl.cv)
misclass.tree(fgl.tr)
misclass.tree(prune.tree(fgl.tr, best=5))

#set.seed(123)
fgl.cv <- cv.tree(fgl.tr,, prune.misclass)
for(i in 2:5)  fgl.cv$dev <- fgl.cv$dev +
    cv.tree(fgl.tr,, prune.misclass)$dev
fgl.cv$dev <- fgl.cv$dev/5
fgl.cv
plot(fgl.cv)
prune.misclass(fgl.tr)

par(mfrow=c(1,1), pty="m")

# End of ch10
