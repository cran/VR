#-*- R -*-

## Script for the Third Edition of `Modern Applied Statistics with S-PLUS'

# Chapter 3   Graphical Output

library(MASS)
postscript(file="ch03.ps", width=8, height=6, pointsize=9)
options(width=65, digits=5)


# 3.2  Basic plotting functions

# data(mdeaths); data(fdeaths)
# library(ts)
lung.deaths <- aggregate(ts.union(mdeaths, fdeaths), 1)
barplot(t(lung.deaths), names=dimnames(lung.deaths)[[1]],
   main="UK deaths from lung disease")
if(interactive())
    legend(locator(1), c("Males", "Females"), fill=c(2,3))
loc <- barplot(t(lung.deaths), names=dimnames(lung.deaths)[[1]],
               angle = c(45, 135), density = 10, col = 1)
total <- apply(lung.deaths, 1, sum)
text(loc, total + par("cxy")[2], total, cex=0.7, xpd=T)

# if(interactive())  brush(as.matrix(hills))

# library(modreg)
# data(topo)
topo.loess <- loess(z ~ x * y, topo, degree=2, span = 0.25)
topo.mar <- list(x = seq(0, 6.5, 0.2), y=seq(0, 6.5, 0.2))
topo.lo <- predict(topo.loess, expand.grid(topo.mar))
topo.lo <- matrix(topo.lo, length(topo.mar$x),length(topo.mar$y))
par(pty="s")       # square plot
contour(topo.mar$x, topo.mar$y, topo.lo, xlab="", ylab="",
   levels = seq(700,1000,25), cex=0.7)
points(topo$x, topo$y)
par(pty="m")
if(F) { # contourplot does not work yet
contourplot(z ~ x * y, con2tr(c(topo.mar, list(z=topo.lo))), aspect=1,
   at = seq(700, 1000, 25), xlab="", ylab="",
   panel = function(x, y, subscripts, ...) {
      panel.contourplot(x, y, subscripts, ...)
      panel.xyplot(topo$x,topo$y, cex=0.5)
   }
)
}


# 3.3  Enhancing plots

# data(wtloss)
attach(wtloss)
oldpar <- par(no.readonly=TRUE)
# alter margin 4; others are default
par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(Days, Weight, type="p",ylab="Weight (kg)")
Wt.lbs <- pretty(range(Weight*2.205))
axis(side=4, at=Wt.lbs/2.205, lab=Wt.lbs, srt=90)
mtext("Weight (lb)", side=4, line=3)
detach()
par(oldpar)


# 3.4  Fine control of graphics

# dataset swiss is done differently in R
data(swiss)
swiss.df <- swiss
attach(swiss.df)
qqnorm(Infant.Mortality)
qqline(Infant.Mortality)

samp <- cbind(Infant.Mortality, matrix(rnorm(47*19), 47, 19))
samp <- apply(scale(samp), 2, sort)
rs <- samp[,1]
xs <- qqnorm(rs, plot=F)$x
env <- t(apply(samp[,-1], 1, range))

matplot(xs, cbind(rs,env), type="pnn", pch=4, mkh=0.06, axes=FALSE,
        xlab="", ylab="")

xyul <- par("usr")
smidge <- min(diff(c(xyul[1], xs, xyul[2])))/2
segments(xs-smidge, env[,1], xs+smidge, env[,1])
segments(xs-smidge, env[,2], xs+smidge, env[,2])

xul <- trunc(10*xyul[1:2])/10
axis(1, at=seq(xul[1], xul[2], by=0.1), labels=FALSE, tck=0.01)
xi <- trunc(xyul[1:2])
axis(1, at=seq(xi[1], xi[2], by=0.5), tck=0.02)
yul <- trunc(5*xyul[3:4])/5
axis(2, at=seq(yul[1], yul[2], by=0.2), labels=FALSE, tck=0.01)
yi <- trunc(xyul[3:4])
axis(2, at=yi[1]:yi[2], tck=0.02)

box(bty="l")          # lower case "L"
# ps.options()$fonts
# R cannot change font family in a plot.
mtext("Quantiles of Standard Normal", side=1, line=2.5, font=3)
mtext("Ri", side=2, line=2, at=yul[2])
detach()
dev.off()


# 3.5  Trellis graphics

library(lattice)
trellis.device(postscript, file="ch03b.ps", width=8, height=6, pointsize=9)

if(F){
# trellis.device()
p1 <- histogram(geyser$waiting)
p1     # plots it on screen

show.settings()
}

# data(hills)
# library(lqs)

xyplot(time ~ dist, data = hills,
  panel = function(x, y, ...) {
     panel.xyplot(x, y, ...)
     panel.lmline(x, y, type="l")
     panel.abline(ltsreg(x, y), lty=3)
#     identify(x, y, row.names(hills)) ## no lattice equivalent
  }
)

# data(michelson)
bwplot(Expt ~ Speed, data=michelson, ylab="Experiment No.",
       main = "Speed of Light Data")
# title("Speed of Light Data") ## fails in lattice

lung.deaths.df <- data.frame(year = rep(1974:1979, 2),
  deaths = c(lung.deaths[, 1], lung.deaths[ ,2]),
  sex = rep(c("Male", "Female"), rep(6,2)))
barchart(year ~ deaths | sex, lung.deaths.df, xlim=c(0, 20000))

splom(~ swiss.df, aspect="fill",
   panel = function(x, y, ...) {
      panel.xyplot(x, y, ...)
      panel.loess(x, y, ...)
   }
)

# data(stormer)
sps <- trellis.par.get("superpose.symbol")
sps$pch <- 1:7
trellis.par.set("superpose.symbol", sps)
xyplot(Time ~ Viscosity, stormer, groups = Wt,
   panel = panel.superpose, type = "b",
   key = list(columns = 3,
       text = list(paste(c("Weight:   ", "", ""),
                         unique(stormer$Wt), "gms")),
       points = Rows(sps, 1:3)
       )
)
rm(sps)

if(F) {
## very slow, incorrect so far
topo.plt <- expand.grid(topo.mar)
topo.plt$pred <- as.vector(predict(topo.loess, topo.plt))
levelplot(pred ~ x * y, topo.plt, aspect=1,
  at = seq(690, 960, 10), xlab="", ylab="",
  panel = function(x, y, subscripts, ...) {
     panel.levelplot(x, y, subscripts, ...)
     panel.xyplot(topo$x,topo$y, cex=0.5, col=1)
  }
)
}

if(F) {
wireframe(pred ~ x * y, topo.plt, aspect=c(1, 0.5), drape=T,
  screen = list(z = -150, x = -60),
  colorkey=list(space="right", height=0.6))
}

# data(crabs)
# library(mva)
lcrabs.pc <- predict(princomp(log(crabs[,4:8])))
crabs.grp <- c("B", "b", "O", "o")[rep(1:4, rep(50,4))]
splom(~ lcrabs.pc[, 1:3], groups = crabs.grp,
   panel = panel.superpose,
   key = list(text = list(c("Blue male", "Blue female",
                           "Orange Male", "Orange female")),
       points = Rows(trellis.par.get("superpose.symbol"), 1:4),
       columns = 4)
  )

sex <- crabs$sex; levels(sex) <- c("Female", "Male")
sp <- crabs$sp; levels(sp) <- c("Blue", "Orange")
splom(~ lcrabs.pc[, 1:3] | sp*sex, cex=0.5, pscales=0)


# data(quine)
Quine <- quine
levels(Quine$Eth) <- c("Aboriginal", "Non-aboriginal")
levels(Quine$Sex) <- c("Female", "Male")
levels(Quine$Age) <- c("primary", "first form",
                       "second form", "third form")
levels(Quine$Lrn) <- c("Average learner", "Slow learner")
bwplot(Age ~ Days | Sex*Lrn*Eth, data=Quine)

bwplot(Age ~ Days | Sex*Lrn*Eth, data=Quine, layout=c(4,2))

stripplot(Age ~ Days | Sex*Lrn*Eth, data=Quine,
   jitter = T, layout = c(4,2))

stripplot(Age ~ Days | Eth*Sex, data=Quine,
   groups = Lrn, jitter=T,
   panel = function(x, y, subscripts, jitter.data=F, ...) {
       if(jitter.data)  y <- jitter(as.numeric(y)) # lattice change
       panel.superpose(x, y, subscripts, ...)
   },
   xlab = "Days of absence",
   between = list(y=1), par.strip.text = list(cex=1.2),
   key = list(columns = 2, text = list(levels(Quine$Lrn)),
       points = Rows(trellis.par.get("superpose.symbol"), 1:2)
       ),
   strip = function(...)
             strip.default(..., strip.names=c(T, T), style=1)
)

# data(fgl)
fgl0 <- fgl[ ,-10] # omit type.
fgl.df <- data.frame(type = rep(fgl$type, 9),
  y = as.vector(as.matrix(fgl0)),
  meas = factor(rep(1:9, rep(214,9)), labels=names(fgl0)))
stripplot(type ~ y | meas, data=fgl.df, scales=list(x="free"),
  strip=function(...) strip.default(style=1, ...), xlab="")

Cath <- equal.count(swiss.df$Catholic, number=2, overlap=0)
Agr5 <- equal.count(swiss.df$Agric, number=5, overlap=0.25)
xyplot(Fertility ~ Education | Agr5 * Cath, data=swiss.df,
          layout=c(2,3), skip = c(F,F,F,F,F,T))

# reorder.factor(Quine$Age, Quine$Days, median)

Cath <- equal.count(swiss.df$Catholic, number=6, overlap=0.25)
xyplot(Fertility ~ Education | Cath, data=swiss.df,
   panel = function(x, y) {
      panel.xyplot(x, y)
      panel.loess(x, y)
   }
)

Cath <- equal.count(swiss.df$Catholic, number=2, overlap=0)
Agr <- equal.count(swiss.df$Agric, number=3, overlap=0.25)
xyplot(Fertility ~ Education | Cath * Agr, data=swiss.df,
   panel = function(x, y) {
      panel.xyplot(x, y)
      panel.loess(x, y)
   }
)

Cath <- equal.count(swiss.df$Cath, number=6, overlap=0.25)
Cath
levels(Cath)
plot(Cath, aspect = 0.3)


# End of ch03
