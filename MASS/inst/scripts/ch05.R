#-*- R -*-

# Chapter 5   Distributions and Data Summaries


# for later use, from section 5.6
perm.t.test <- function(d) {
# ttest is function(x) mean(x)/sqrt(var(x)/length(x))
    binary.v <- function(x, digits) {
       if(missing(digits)) {
           mx <- max(x)
           digits <- if(mx > 0) 1 + floor(log(mx, base = 2)) else 1
       }
       ans <- 0:(digits - 1)
       lx <- length(x)
       x <- rep(x, rep(digits, lx))
       x <- (x %/% 2^ans) %% 2
       dim(x) <- c(digits, lx)
       x
    }
    digits <- length(d)
    n <- 2^digits
    x <- d * 2 * (binary.v(1:n, digits) - 0.5)
    mx <- matrix(1/digits, 1, digits) %*% x
    s <- matrix(1/(digits - 1), 1, digits)
    vx <- s %*% (x - matrix(mx, digits, n, byrow=T))^2
    as.vector(mx/sqrt(vx/digits))
}
library(MASS)
options(width=65, digits=5, height=9999)
postscript(file="ch05.ps", width=8, height=6, pointsize=9)

rm(A, B) # precautionary clear-out
data(shoes)
attach(shoes)
tperm <- perm.t.test(B - A) # see section 5.6
detach()


# 5.1  Probability distributions

x <- rt(250, 9)
qqnorm(x); qqline(x)

library(lattice)
trellis.device(postscript, file="ch05b.ps", width=8, height=6, pointsize=9)
qqmath(~ x, distribution=qnorm, aspect="xy",
   prepanel = prepanel.qqmathline,
   panel = function(x, y, ...) {
      panel.qqmathline(y, distribution=qnorm, ...)
      panel.qqmath(x, y, ...)
   },
   xlab = "Quantiles of Standard Normal"
)
dev.off()
grid.stop()

# 5.2  Generating random data

contam <- rnorm( 100, 0, (1 + 2*rbinom(100, 1, 0.05)) )


# 5.3  Data summaries

data(geyser)
data(chem)
data(abbey)
par(mfrow=c(2,3))
hist.scott(geyser$duration, xlab="duration")
hist.scott(chem)
hist.scott(tperm)
hist.FD(geyser$duration, xlab="duration")
hist.FD(chem)
hist.FD(tperm)
par(mfrow=c(1,1))

data(swiss)
swiss.fertility <- swiss[, 1]
stem(swiss.fertility)
stem(chem)
stem(abbey)
stem(abbey, scale=0.4) # different in R

par(mfrow=c(1,2))
boxplot(chem, sub="chem", range=0.5)
boxplot(abbey, sub="abbey")
par(mfrow=c(1,1))
# fgl.df is from Chapter 3
data(fgl)
fgl0 <- fgl[ ,-10] # omit type.
fgl.df <- data.frame(type = rep(fgl$type, 9),
   y = as.vector(as.matrix(fgl0)),
   meas = factor(rep(1:9, rep(214,9)), labels=names(fgl0)))
#bwplot(type ~ y | meas, data=fgl.df, scales=list(x="free"),
#   strip=function(...) strip.default(..., style=1), xlab="")


# 5.4  Classical univariate statistics

attach(shoes)
t.test(A, mu = 10)
t.test(A)$conf.int
wilcox.test(A, mu = 10)
var.test(A, B)
t.test(A, B)
t.test(A, B, var.equal=F)
wilcox.test(A, B)
t.test(A, B, paired=T)
wilcox.test(A, B, paired=T)
detach()

par(mfrow=c(1,2))
truehist(tperm, xlab="diff")
x <- seq(-4,4, 0.1)
lines(x, dt(x,9))
#cdf.compare(tperm, distribution="t", df=9)
sres <- c(sort(tperm), 4)
yres <- (0:1024)/1024
plot(sres, yres, type="S", xlab="diff", ylab="")
lines(x, pt(x,9), lty=3)
legend(-5, 1.05, c("Permutation dsn","t_9 cdf"), lty=c(1,3))
par(mfrow=c(1,1))


# 5.5  Robust summaries

sort(chem)
mean(chem)
median(chem)
#location.m(chem)
#location.m(chem, psi.fun="huber")
mad(chem)
#scale.tau(chem)
#scale.tau(chem, center=3.68)
unlist(huber(chem))
unlist(hubers(chem))

sort(abbey)
mean(abbey)
median(abbey)
#location.m(abbey)
#location.m(abbey, psi.fun="huber")
unlist(hubers(abbey))
unlist(hubers(abbey, k=2))
unlist(hubers(abbey, k=1))


# 5.6  Density estimation

attach(geyser)
par(mfrow=c(2,3))
truehist(duration, h=0.5, x0=0.0, xlim=c(0, 6), ymax=0.7)
truehist(duration, h=0.5, x0=0.1, xlim=c(0, 6), ymax=0.7)
truehist(duration, h=0.5, x0=0.2, xlim=c(0, 6), ymax=0.7)
truehist(duration, h=0.5, x0=0.3, xlim=c(0, 6), ymax=0.7)
truehist(duration, h=0.5, x0=0.4, xlim=c(0, 6), ymax=0.7)

breaks <- seq(0, 5.9, 0.1)
counts <- numeric(length(breaks))
for(i in (0:4)) counts[i+(1:55)] <- counts[i+(1:55)] +
    rep(hist(duration, breaks=0.1*i + seq(0, 5.5, 0.5),
    prob=T, plot=F)$intensities, rep(5,11))
plot(breaks+0.05, counts/5, type="l", xlab="duration",
    ylab="averaged", bty="n", xlim=c(0, 6), ylim=c(0, 0.7))

par(mfrow=c(2,2))
x <- seq(-5, 5, 0.1)
plot(c(-5, 5), c(0, 0.45), type="n", bty="l",
  sub="default", xlab="", ylab="")
lines(x, dt(x,9), lty=2); rug(jitter(tperm))
lines(density(tperm, n=200, from=-5, to=5))
plot(c(-5, 5), c(0, 0.45), type="n", bty="l",
  sub="width=0.2", xlab="", ylab="")
lines(x, dt(x,9), lty=2); rug(jitter(tperm))
lines(density(tperm, n=200, width=0.2, from=-5, to=5))
plot(c(-5, 5), c(0, 0.45), type="n", bty="l",
  sub="width=0.5", xlab="", ylab="")
lines(x, dt(x,9), lty=2); rug(jitter(tperm))
lines(density(tperm, n=200, width=0.5, from=-5, to=5))
plot(c(-5, 5), c(0, 0.45), type="n", bty="l",
  sub="width=1.5", xlab="", ylab="")
lines(x, dt(x,9), lty=2); rug(jitter(tperm))
lines(density(tperm, n=200, width=1.5, from=-5, to=5))

par(mfrow=c(2,2))
truehist(duration, nbins=15, xlim=c(0.5,6), ymax=1.2)
lines(density(duration, n=200))
bandwidth.nrd(duration)
lines(density(duration, width=1.5565, n=200), lty=3)
bandwidth.nrd(tperm)

c(width.SJ(duration, method="dpi"), width.SJ(duration))
truehist(duration, nbins=15, xlim=c(0.5,6), ymax=1.2)
lines(density(duration, width=0.57, n=200), lty=1)
lines(density(duration, width=0.36, n=200), lty=3)
c( ucv(tperm), bcv(tperm), width.SJ(tperm) )
par(mfrow=c(1,1))

data(galaxies)
gal <- galaxies/1000
c(width.SJ(gal, method="dpi"), width.SJ(gal))
plot(x=c(0, 40), y=c(0, 0.3), type="n", bty="l",
     xlab="velocity of galaxy (1000km/s)", ylab="density")
rug(gal)
lines(density(gal, width=3.25, n=200), lty=1)
lines(density(gal, width=2.56, n=200), lty=3)
median(gal)

library(logspline)
x <- seq(5, 40, length=500)
lines(x, dlogspline(x, logspline.fit(gal)), lty=2)



plot(duration, waiting, xlim=c(0.5,6), ylim=c(40,100))
f1 <- kde2d(duration, waiting, n=50, lims=c(0.5,6,40,100))
image(f1, zlim = c(0, 0.05), col=grey(128:0/128))
# levelplot(z ~ x*y, con2tr(f1),
#    at = seq(0, 0.07, 0.001), colorkey=F,
#    col.regions = rev(trellis.par.get("regions")$col))
f2 <- kde2d(duration, waiting, n=50, lims=c(0.5,6,40,100),
   h = c(width.SJ(duration), width.SJ(waiting)) )
#levelplot(z ~ x*y, con2tr(f2),
#   xlab="duration", ylab="waiting",
#   at = seq(0, 0.07, 0.001), colorkey=F,
#   col.regions = rev(trellis.par.get("regions")$col))
#wireframe(z ~ x*y, con2tr(f2),
#   aspect = c(1, 0.5), screen=list(z=20, x=-60), zoom=1.2)
image(f2, zlim = c(0, 0.05), col=grey(128:0/128))
persp(f2, phi=30, theta=20, d=5)

plot(duration[-272], duration[-1], xlim=c(0.5, 6),
    ylim=c(1, 6),xlab="previous duration", ylab="duration")
f1 <- kde2d(duration[-272], duration[-1],
   h=rep(1.5, 2), n=50, lims=c(0.5,6,0.5,6))
contour(f1 ,xlab="previous duration",
    ylab="duration", levels = c(0.05, 0.1, 0.2, 0.4) )
f1 <- kde2d(duration[-272], duration[-1],
   h=rep(0.6, 2), n=50, lims=c(0.5,6,0.5,6))
contour(f1 ,xlab="previous duration",
    ylab="duration", levels = c(0.05, 0.1, 0.2, 0.4) )
f1 <- kde2d(duration[-272], duration[-1],
   h=rep(0.4, 2), n=50, lims=c(0.5,6,0.5,6))
contour(f1 ,xlab="previous duration",
    ylab="duration", levels = c(0.05, 0.1, 0.2, 0.4) )
detach("geyser")


# 5.7  Bootstrap and permutation methods

density(gal, n=1, from=20.833, to=20.834, width=3.2562)$y
density(gal, n=1, from=20.833, to=20.834, width=2.5655)$y
1/(2*sqrt(length(gal))*0.13)
set.seed(101); m <- 1000
res <- numeric(m)
for (i in 1:m) res[i] <- median(sample(gal, replace=T))
mean(res - median(gal))
sqrt(var(res))

truehist(res, h=0.1)
lines(density(res, width=width.SJ(res, method="dpi"), n=200))
quantile(res, p=c(0.025, 0.975))
x <- seq(19.5, 22.5, length=500)
lines(x, dlogspline(x, logspline.fit(res)), lty=3)

library(boot)
set.seed(101)
gal.boot <- boot(gal, function(x,i) median(x[i]), R=1000)
gal.boot
plot(gal.boot)

if(F){
if(version$major >= 4) {
  gal.bt <- bootstrap(gal, median, seed=101, B=1000)
  print(summary(gal.bt))
  print(limits.emp(gal.bt))
  print(limits.bca(gal.bt))
  plot(gal.bt)
  qqnorm(gal.bt)
}
}

sim.gen  <- function(data, mle) {
  n <- length(data)
  data[sample(n, replace=T)]  + mle*rnorm(n)
}
gal.boot2 <- boot(gal, median, R=1000,
  sim="parametric", ran.gen=sim.gen, mle=0.5)
boot.ci(gal.boot2, conf=c(0.90, 0.95),
        type=c("norm","basic","perc"))

attach(shoes)
t.test(B - A)
shoes.boot <- boot(B-A, function(x,i) mean(x[i]), R=1000)
boot.ci(shoes.boot, type = c("norm", "basic", "perc", "bca"))
mean.fun <- function(d, i) {
  n <- length(i)
  c(mean(d[i]), (n-1)*var(d[i])/n^2)
}
shoes.boot2 <- boot(B-A, mean.fun, R=1000)
boot.ci(shoes.boot2, type = "stud")

d <- B - A
ttest <- function(x) mean(x)/sqrt(var(x)/length(x))
n <- 1000
res <- numeric(n)
for(i in 1:n) res[i] <- ttest(x <- d*sign(runif(10)-0.5))


# End of ch05
