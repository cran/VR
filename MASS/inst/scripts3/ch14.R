#-*- R -*-

## Script for the Third Edition of `Modern Applied Statistics with S-PLUS'

# Chapter 14   Spatial Statistics

postscript(file="ch14.ps", width=8, height=8, pointsize=9)
options(width=65, digits=5)

library(MASS)
library(spatial)

# 14.1  Spatial interpolation and smoothing

data(topo)
par(mfrow=c(2,2), pty="s")
topo.ls <- surf.ls(2, topo)
trsurf <- trmat(topo.ls, 0, 6.5, 0, 6.5, 30)
contour(trsurf, levels=seq(600,1000,25), xlab="", ylab="")
points(topo)
title("Degree=2")
topo.ls <- surf.ls(3, topo)
trsurf <- trmat(topo.ls, 0, 6.5, 0, 6.5, 30)
contour(trsurf, levels=seq(600,1000,25), xlab="", ylab="")
points(topo)
title("Degree=3")
topo.ls <- surf.ls(4, topo)
trsurf <- trmat(topo.ls, 0, 6.5, 0, 6.5, 30)
contour(trsurf, levels=seq(600,1000,25), xlab="", ylab="")
points(topo)
title("Degree=4")
topo.ls <- surf.ls(6, topo)
trsurf <- trmat(topo.ls, 0, 6.5, 0, 6.5, 30)
contour(trsurf, levels=seq(600,1000,25), xlab="", ylab="")
points(topo)
title("Degree=6")

topo.ls <- surf.ls(4, topo)
trsurf <- trmat(topo.ls, 0, 6.5, 0, 6.5, 30)
# trsurf[c("x", "y")] <- expand.grid(x=trsurf$x, y=trsurf$y)
# plt1 <- levelplot(z ~ x * y, trsurf, aspect=1,
#            at = seq(650, 1000, 10),  xlab = "", ylab = "")
# plt2 <- wireframe(z ~ x * y, trsurf, aspect=c(1, 0.5),
#            screen = list(z = -30, x = -60))
# print(plt1, position = c(0, 0, 0.5, 1), more=T)
# print(plt2, position = c(0.45, 0, 1, 1))
image(trsurf, col=grey(0:128/128))
points(topo)

library(modreg)
dif <- 0.1
par(mfcol=c(2,2), pty="s")
topo.loess <- loess(z ~ x * y, topo, degree=2, span = 0.25, normalize=F)
topo.mar <- list(x = seq(0, 6.5, dif), y=seq(0, 6.5, dif))
topo.lo <- predict(topo.loess, expand.grid(topo.mar), se=T)
topo.lo$fit <- matrix(topo.lo$fit, length(topo.mar$x), length(topo.mar$y))
contour(topo.mar$x,topo.mar$y,topo.lo$fit, levels = seq(700,1000,25),
        xlab="fit", ylab="")
points(topo)
topo.lo$se.fit <- matrix(topo.lo$se.fit, length(topo.mar$x), length(topo.mar$y))
contour(topo.mar$x,topo.mar$y,topo.lo$se.fit, levels = seq(5, 25, 5),
        xlab="standard error", ylab="")
title("Loess degree = 2")
points(topo)

topo.loess <- loess(z ~ x * y, topo, degree=1, span = 0.25, normalize=F,
                    xlab="", ylab="")
topo.lo <- predict(topo.loess, expand.grid(topo.mar), se=T)
topo.lo$fit <- matrix(topo.lo$fit, length(topo.mar$x), length(topo.mar$y))
contour(topo.mar$x,topo.mar$y,topo.lo$fit, levels = seq(700,1000,25),
        xlab="fit", ylab="")
points(topo)
topo.lo$se.fit <- matrix(topo.lo$se.fit, length(topo.mar$x), length(topo.mar$y))
contour(topo.mar$x,topo.mar$y,topo.lo$se.fit, levels = seq(5, 25, 5),
        xlab="standard error", ylab="")
title("Loess degree = 1")
points(topo)

library(akima)
par(mfrow=c(1,2), pty="s")
contour(interp.old(topo$x, topo$y, topo$z),xlab="interp default",
        ylab="", levels = seq(600,1000,25))
points(topo)
topo.mar <- list(x = seq(0, 6.5, 0.1), y=seq(0, 6.5, 0.1))
contour(interp.old(topo$x, topo$y, topo$z, topo.mar$x, topo.mar$y,
   ncp=4, extrap=T), xlab="interp", ylab="",
   levels = seq(600,1000,25))
points(topo)



# 14.2  Kriging

par(mfrow=c(2,2), pty="s")
topo.ls <- surf.ls(2, topo)
trsurf <- trmat(topo.ls, 0, 6.5, 0, 6.5, 30)
contour(trsurf, levels=seq(600, 1000, 25), xlab="", ylab="")
points(topo)
title("LS trend surface")
topo.gls <- surf.gls(2, expcov, topo, d=0.7)
trsurf <- trmat(topo.gls, 0, 6.5, 0, 6.5, 30)
contour(trsurf, levels=seq(600, 1000, 25), xlab="", ylab="")
points(topo)
title("GLS trend surface")

prsurf <- prmat(topo.gls, 0, 6.5, 0, 6.5, 50)
contour(prsurf, levels=seq(600, 1000, 25), xlab="", ylab="")
points(topo)
title("Kriging prediction")
sesurf <- semat(topo.gls, 0, 6.5, 0, 6.5, 30)
contour(sesurf, levels=c(20, 25), xlab="", ylab="")
points(topo)
title("Kriging s.e.")

par(mfrow=c(2,2), pty="m")
topo.kr <- surf.ls(2, topo)
correlogram(topo.kr, 25)
d <- seq(0, 7, 0.1)
lines(d, expcov(d, 0.7))
variogram(topo.kr, 25)

topo.kr <- surf.gls(2, expcov, topo, d=0.7)
correlogram(topo.kr, 25)
lines(d, expcov(d, 0.7))
lines(d, gaucov(d, 1.0, 0.3), lty=3) # try nugget effect

topo.kr <- surf.ls(0, topo)
correlogram(topo.kr, 25)
lines(d, gaucov(d, 2, 0.05))

par(mfrow=c(2,2), pty="s")
topo.kr <- surf.gls(2, gaucov, topo, d=1, alph=0.3)
prsurf <- prmat(topo.kr, 0, 6.5, 0, 6.5, 50)
contour(prsurf, levels=seq(600, 1000, 25), xlab="fit", ylab="")
points(topo)
sesurf <- semat(topo.kr, 0, 6.5, 0, 6.5, 25)
contour(sesurf, levels=c(15, 20, 25), xlab="standard error", ylab="")
points(topo)

topo.kr <- surf.gls(0, gaucov, topo, d=2, alph=0.05, nx=10000)
prsurf <- prmat(topo.kr, 0, 6.5, 0, 6.5, 50)
contour(prsurf, levels=seq(600, 1000, 25), xlab="fit", ylab="")
points(topo)
sesurf <- semat(topo.kr, 0, 6.5, 0, 6.5, 25)
contour(sesurf, levels=c(15, 20, 25), xlab="standard error", ylab="")
points(topo)



# 14.3  Point process analysis

library(spatial)
pines <- ppinit("pines.dat")
par(mfrow=c(2,2), pty="s")
plot(pines, xlim=c(0,10), ylim=c(0,10), xlab="", ylab="",
     xaxs="i", yaxs="i")
plot(Kfn(pines,5), type="s", xlab="distance", ylab="L(t)")
lims <- Kenvl(5, 100, Psim(72))
lines(lims$x, lims$l, lty=2)
lines(lims$x, lims$u, lty=2)
ppregion(pines)
plot(Kfn(pines,1.5), type="s", xlab="distance", ylab="L(t)")
lims <- Kenvl(1.5, 100, Strauss(72, 0.2, 0.7))
lines(lims$x, lims$a, lty=2)
lines(lims$x, lims$l, lty=2)
lines(lims$x, lims$u, lty=2)
pplik(pines, 0.7)
lines(Kaver(1.5, 100, Strauss(72, 0.15, 0.7)), lty=3)

# End of ch14
