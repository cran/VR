#-*- R -*-

## Script for the Third Edition of `Modern Applied Statistics with S-PLUS'

# Chapter 8   Non-linear Models

library(MASS)
library(nls)
postscript(file="ch08.ps", width=8, height=6, pointsize=9)
options(width=65, digits=5)

# 8.1 An introductory example

data(wtloss)
attach(wtloss)
# alter margin 4; others are default
oldpar <- par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(Days, Weight, type="p", ylab="Weight (kg)")
Wt.lbs <- pretty(range(Weight*2.205))
axis(side=4, at=Wt.lbs/2.205, lab=Wt.lbs, srt=90)
mtext("Weight (lb)", side=4, line=3)
par(oldpar) # restore settings
detach()


# 8.2  Fitting non-linear regression models

wtloss.st <- c(b0=90, b1=95, th=120)
wtloss.fm <- nls(Weight ~ b0 + b1*2^(-Days/th),
   data = wtloss, start = wtloss.st, trace = T)
wtloss.fm

expn <- function(b0, b1, th, x)
{
    temp <- 2^(-x/th)
    model.func <- b0 + b1 * temp
    Z <- cbind(1, temp, (b1 * x * temp * log(2))/th^2)
    dimnames(Z) <- list(NULL, c("b0","b1","th"))
    attr(model.func, "gradient") <- Z
    model.func
}

wtloss.gr <- nls(Weight ~ expn(b0, b1, th, Days),
   data = wtloss, start = wtloss.st, trace = T)

## R needs a different syntax here
expn1 <- deriv(y ~ b0 + b1 * 2^(-x/th), c("b0", "b1", "th"),
               c("b0", "b1", "th", "x"))
expn1

negexp <- selfStart(model = ~ b0 + b1*exp(-x/th),
    initial = negexp.SSival, parameters = c("b0", "b1", "th"),
    template = function(x, b0, b1, th) {})

wtloss.ss <- nls(Weight ~ negexp(Days, B0, B1, theta),
                 data=wtloss, trace = T)


# 8.3  Non-linear fitted model objects and method functions

summary(wtloss.gr)
deviance(wtloss.gr)
vcov(wtloss.gr)

data(muscle)
A <- model.matrix(~ Strip - 1, data=muscle)
rats.nls1 <- nls(log(Length) ~ cbind(A, rho^Conc),
   data = muscle, start = c(rho=0.1), algorithm="plinear")
B <- coef(rats.nls1)
B

st <- list(alpha = B[2:22], beta = B[23], rho = B[1])
rats.nls2 <- nls(log(Length) ~ alpha[Strip] + beta*rho^Conc,
                  data = muscle, start = st)

attach(muscle)
Muscle <- expand.grid(Conc = sort(unique(Conc)),
                     Strip = levels(Strip))
Muscle$Yhat <- predict(rats.nls2, Muscle)
Muscle$logLength <- as.numeric(rep(NA, nrow(Muscle)))
ind <- match(paste(Strip, Conc),
            paste(Muscle$Strip, Muscle$Conc))
Muscle$logLength[ind] <- log(Length)
detach()

if(F) { # no Trellis
xyplot(Yhat ~ Conc | Strip, Muscle, as.table = T,
  ylim = range(c(Muscle$Yhat, Muscle$logLength), na.rm=T),
  subscripts = T, xlab = "Calcium Chloride concentration (mM)",
  ylab = "log(Length in mm)", panel =
  function(x, y, subscripts, ...) {
     lines(spline(x, y))
     panel.xyplot(x, Muscle$logLength[subscripts], ...)
  })
}

if(F) { ## prior to 1.1.0
coplot(seq(0.8, 4, len=126) ~ Conc | Strip, Muscle, show.given=FALSE,
       xlab = "Calcium Chloride concentration (mM)",
       ylab = "log(Length in mm)", panel = function(x, y, ...) {
           ind <- round(1+125*(y-0.8)/3.2)
           lines(spline(x, Muscle$Yhat[ind]))
           points(x, Muscle$logLength[ind])
       })
}

## from 1.1.0
coplot(Yhat ~ Conc | Strip, Muscle, show.given=FALSE,
       xlab = "Calcium Chloride concentration (mM)",
       ylab = "log(Length in mm)", subscripts = TRUE,
       panel = function(x, y, subscripts, ...) {
           points(x, y)
           lines(spline(x, Muscle$Yhat[subscripts]))
       })

# 8.5  Confidence intervals for parameters

expn2 <- deriv(~b0 + b1*((w0 - b0)/b1)^(x/d0),
        c("b0","b1","d0"), function(b0, b1, d0, x, w0) {})

wtloss.init <- function(obj, w0) {
  p <- coef(obj)
  d0 <-  - log((w0 - p["b0"])/p["b1"])/log(2) * p["th"]
  c(p[c("b0", "b1")], d0 = as.vector(d0))
}

out <- NULL
w0s <- c(110, 100, 90)
for(w0 in w0s) {
    fm <- nls(Weight ~ expn2(b0, b1, d0, Days, w0),
              wtloss, start = wtloss.init(wtloss.gr, w0))
    out <- rbind(out, c(coef(fm)["d0"], confint(fm, "d0")))
  }
dimnames(out)[[1]] <- paste(w0s,"kg:")
out

data(stormer)
attach(stormer)
fm0 <- lm(Wt*Time ~ Viscosity + Time - 1,  data=stormer)
b0 <- coef(fm0)
names(b0) <- c("b1", "b2")
b0
storm.fm <- nls(Time ~ b1*Viscosity/(Wt-b2), data=stormer, start=b0,
          trace=T)

bc <- coef(storm.fm)
se <- sqrt(diag(vcov(storm.fm)))
dv <- deviance(storm.fm)

par(pty = "s")
b1 <- bc[1] + seq(-3*se[1], 3*se[1], length = 51)
b2 <- bc[2] + seq(-3*se[2], 3*se[2], length = 51)
bv <- expand.grid(b1, b2)
ssq <- function(b)
      sum((Time - b[1] * Viscosity/(Wt-b[2]))^2)
dbetas <- apply(bv, 1, ssq)
cc <- matrix(Time - rep(bv[,1],rep(23, 2601)) *
      Viscosity/(Wt - rep(bv[,2], rep(23, 2601))), 23)
dbetas <- matrix(drop(rep(1, 23) %*% cc^2), 51)
fstat <- matrix( ((dbetas - dv)/2) / (dv/21), 51, 51)
qf(0.95, 2, 21)
plot(b1, b2, type="n")
lev <- c(1, 2, 5, 7, 10, 15, 20)
contour(b1, b2, fstat, levels=lev, labex=0.75, lty=2, add=T)
contour(b1, b2, fstat, levels=qf(0.95,2,21), add=T, labex=0)
text(31.6, 0.3, labels="95% CR", adj=0, cex=0.75)
points(bc[1], bc[2], pch=3, mkh=0.1)
detach()

library(boot)
storm.fm <- nls(Time ~ b*Viscosity/(Wt - c), stormer,
                start = c(b=29.401, c=2.2183))
summary(storm.fm)$parameter
st <- cbind(stormer, fit=fitted(storm.fm))
storm.bf <- function(rs, i) {
#  st <- st                              # for S-PLUS 5.x
  st$Time <-  st$fit + rs[i]
  coef(nls(Time ~ (b * Viscosity)/(Wt - c), st,
           start = coef(storm.fm)))
}
rs <- scale(resid(storm.fm), scale = F) # remove the mean
storm.boot <- boot(rs, storm.bf, R = 4999) # pretty slow
storm.boot
boot.ci(storm.boot, index=1,
        type=c("norm", "basic", "perc", "bca"))
boot.ci(storm.boot, index=2,
        type=c("norm", "basic", "perc", "bca"))


# 8.5  Assessing the linear approximation

opar <- par(pty="m", mfrow=c(1,3))
plot(profile(update(wtloss.gr, trace=F)))
par(opar)


# 8.6  Constrained non-linear regression

data(whiteside)
attach(whiteside)
Gas <- Gas[Insul=="Before"]
Temp <- -Temp[Insul=="Before"]
#nnls.fit(cbind(1, -1, Temp), Gas)
# can use box-constrained optimizer
fn <- function(par) sum((Gas - par[1] - par[2]*Temp)^2)
optim(rep(0,2), fn, lower=0, method="L-BFGS-B")$par
rm(Gas, Temp)
detach()

data(wtloss)
attach(wtloss)
optim(c(90,95,120),
      function(x) sum((Weight-x[1]-x[2]*2^(-Days/x[3]))^2),
      lower=rep(0,3), method="L-BFGS-B")$par
detach()

if(F) { # have none of these
wtloss.r <- function(x, Weight, Days)
    Weight - x[1] - x[2] * 2^(-Days/x[3])
wtloss.rg <- function(x, Weight, Days)
{
    temp <- 2^(-Days/x[3])
    -cbind(1, temp, x[2]*Days*temp*log(2)/x[3]^2)
}
wtloss.nl <- nlregb(nrow(wtloss), c(90, 95, 120),
   wtloss.r,  wtloss.rg, lower = rep(0,3),
   Weight = wtloss$Weight, Days = wtloss$Days)

vcov1 <- function(object)
{
   gr <- object$jacobian
   df <- length(object$resid) - length(object$param)
   sum(object$resid^2)/df * solve(t(gr) %*% gr)
}

sqrt(diag(vcov1(wtloss.nl)))
sqrt(diag(vcov.nlregb(wtloss.nl, method="Fisher")))
sqrt(diag(vcov.nlregb(wtloss.nl, method="observed")))
sqrt(diag(vcov.nlregb(wtloss.nl, method="Huber")))

wtloss.nl0 <-  nlregb(nrow(wtloss), c(90,95,120),
    wtloss.r, lower = rep(0,3),
    Weight = wtloss$Weight, Days = wtloss$Days)
sqrt(diag(vcov.nlregb(wtloss.nl0)))
}


# 8.7  General optimization and maximum likelihood estimation

data(geyser)
attach(geyser)
truehist(waiting, xlim=c(35,110), ymax=0.04, h=5)
width.SJ(waiting)
wait.dns <- density(waiting, n=200, width=10.24)
lines(wait.dns, lty=2)

if(F) { # have none of these
lmix2 <- deriv3(
     ~ -log(p*dnorm((x-u1)/s1)/s1 + (1-p)*dnorm((x-u2)/s2)/s2),
     c("p", "u1", "s1", "u2", "s2"),
     function(x, p, u1, s1, u2, s2) NULL)

p0 <- c(p=mean(waiting < 70), u1=50, s1=5, u2=80, s2=5)
p0
tr.ms <- function(info, theta, grad, scale, flags, fit.pars)
{
    cat(round(info[3], 3), ":", signif(theta), "\n")
    invisible()
}

wait.mix2 <- ms(~ lmix2(waiting, p, u1, s1, u2, s2),
   start=p0, data = geyser, trace = tr.ms)

dmix2 <- function(x, p, u1, s1, u2, s2)
            p * dnorm(x, u1, s1) + (1-p) * dnorm(x, u2, s2)
cf <- coef(wait.mix2)
attach(structure(as.list(cf), names = names(cf)))
wait.fdns <- list(x = wait.dns$x,
                  y = dmix2(wait.dns$x, p, u1, s1, u2, s2))
lines(wait.fdns)
par(usr = c(0,1,0,1))
legend(0.1, 0.9, c("Normal mixture", "Nonparametric"),
    lty = c(1,2), bty = "n")


pmix2 <- deriv(~ p*pnorm((x-u1)/s1) + (1-p)*pnorm((x-u2)/s2),
               "x", function(x, p, u1, s1, u2, s2) {})
pr0 <- (seq(along = waiting) - 0.5)/length(waiting)
x0 <- x1 <- as.vector(sort(waiting)) ; del <- 1; i <- 0
while((i <- 1 + 1) < 10 && abs(del) > 0.0005) {
    pr <- pmix2(x0, p, u1, s1, u2, s2)
    del <- (pr - pr0)/attr(pr, "gradient")
    x0 <- x0 - 0.5*del
    cat(format(del <- max(abs(del))), "\n")
}
detach()

par(pty = "s")
plot(x0, x1, xlim = range(x0, x1), ylim = range(x0, x1),
    xlab = "Model quantiles", ylab = "Waiting time")
abline(0,1)
par(pty = "m")

vmat <- summary(wait.mix2)$Information
cbind(coef(wait.mix2), sqrt(diag(vmat)))


lmix2r <- deriv3(
     ~ -log((exp(a+b*y)*dnorm((x-u1)/s1)/s1 +
             dnorm((x-u2)/s2)/s2) / (1+exp(a+b*y)) ),
     c("a", "b", "u1", "s1", "u2", "s2"),
     function(x, y, a, b, u1, s1, u2, s2) NULL)

p1 <- wait.mix2$par
tmp <- as.vector(p1[1])
p2 <- c(a=log(tmp/(1-tmp)), b=0, p1[-1])

wait.mix2r <- ms(~ lmix2r(waiting[-1], duration[-299], a, b, u1, s1, u2, s2),
                 start = p2, data = geyser, trace = tr.ms)

grid <- expand.grid(x=seq(1.5, 5.5, 0.1), y=seq(40, 110, 0.5))
grid$z <- exp(-lmix2r(grid$y, grid$x, 16.14, -5.74, 55.14, 5.663, 81.09, 6.838))

levelplot(z ~ x*y, grid, colorkey=F, at = seq(0, 0.075, 0.001),
          panel= function(...) {
            panel.levelplot(...)
            points(duration[-299], waiting[-1])
          }, xlab="previous duration", ylab="wait",
   col.regions = rev(trellis.par.get("regions")$col))
}

mix.f <- function(p)
{
   e <- p[1]*dnorm((waiting-p[2])/p[3])/p[3] +
        (1-p[1])*dnorm((waiting-p[4])/p[5])/p[5]
   if(any(e <= 0)) Inf else -sum(log(e))
}
waiting.init <- c(mean(waiting < 70), 50, 5, 80, 5)
nlm(mix.f, waiting.init, print.level=1)

mix.obj <- function(p, x)
{
   e <- p[1]*dnorm((x-p[2])/p[3])/p[3] +
        (1-p[1])*dnorm((x-p[4])/p[5])/p[5]
   if(any(e <= 0)) Inf else -sum(log(e))
}
optim(waiting.init, mix.obj, x = waiting)
optim(waiting.init, mix.obj, method="BFGS", x = waiting)

mix.nl0 <- optim(waiting.init, mix.obj, method="L-BFGS-B",
                 lower = c(0, -Inf, 0, -Inf, 0),
                 upper = c(1, rep(Inf, 4)), x = waiting)

# mix.nl0 <- nlminb(waiting.init,  mix.obj,
#    scale = c(10, rep(1,4)), lower = c(0, -Inf, 0, -Inf, 0),
#    upper = c(1, rep(Inf, 4)), x = waiting)

if(F) {
lmix2a <- deriv(
     ~ -log(p*dnorm((x-u1)/s1)/s1 + (1-p)*dnorm((x-u2)/s2)/s2),
     c("p", "u1", "s1", "u2", "s2"),
     function(x, p, u1, s1, u2, s2) NULL)

mix.gr <- function(p, x)
{
   u1 <- p[2]; s1 <- p[3]; u2 <- p[4]; s2 <- p[5]; p <- p[1]
   e <- lmix2a(x, p, u1, s1, u2, s2)
   rep(1, length(x)) %*% attr(e, "gradient")
}
mix.nl1 <- nlminb(waiting.init, mix.obj, mix.gr,
   scale = c(10, rep(1,4)), lower = c(0, -Inf, 0, -Inf, 0),
   upper = c(1, rep(Inf, 4)), x = waiting)

mix.grh <- function(p, x)
{
   e <- lmix2(x, p[1], p[2], p[3], p[4], p[5])
   g <- attr(e, "gradient")
   g <- rep(1, length(x)) %*% g
   H <- apply(attr(e, "hessian"), c(2,3), sum)
   list(gradient=g, hessian=H[row(H) <= col(H)])
}
mix.nl2 <- nlminb(waiting.init, mix.obj, mix.grh, T,
   scale = c(10, rep(1,4)), lower = c(0, -Inf, 0, -Inf, 0),
   upper = c(1, rep(Inf, 4)), x = waiting)

sqrt(diag(vcov.nlminb(mix.nl0)))
sqrt(diag(vcov.nlminb(mix.nl1)))
sqrt(diag(vcov.nlminb(mix.nl2)))
}
detach(geyser)

AIDSfit <- function(y, z, start=rep(mean(y), ncol(z)), ...)
{
  deviance <- function(beta, y, z) {
      mu <- z %*% beta
      2 * sum(mu - y - y*log(mu/y)) }
  grad <- function(beta, y, z) {
      mu <- z %*% beta
      2 * t(1 - y/mu) %*% z }
  optim(start, deviance, grad, lower = 0, y = y, z = z,
        method="L-BFGS-B", ...)
}

Y <- scan(n=13)
12 14 33 50 67 74 123 141 165 204 253 246 240

library(nnet) # for class.ind
s <- seq(0, 13.999, 0.01); tint <- 1:14
X <- expand.grid(s, tint)
Z <- matrix(pweibull(pmax(X[,2] - X[,1],0), 2.5, 10),length(s))
Z <- Z[,2:14] - Z[,1:13]
Z <- t(Z) %*% class.ind(factor(floor(s/2))) * 0.01
round(AIDSfit(Y, Z)$par)
rm(s, X, Y, Z)


# 8.8 Non-linear mixed effects models

library(nlme)
options(contrasts = c("contr.treatment", "contr.poly"))
data(Sitka)
sitka.nlme <- nlme(size ~ A + B * (1 - exp(-(Time-100)/C)),
   fixed = list(A ~ treat, B ~ treat, C ~ 1),
   random =   A + B ~ 1 | tree, data = Sitka,
   start  = list(fixed = c(2, 0, 4, 0, 80)),
   method = "ML", verbose = T)

summary(sitka.nlme)

sitka.nlme2 <- update(sitka.nlme,
    fixed = list(A ~ 1, B ~ 1, C ~ 1),
    start = list(fixed=c(2.4, 3.6, 82)))
summary(sitka.nlme2)
anova(sitka.nlme2, sitka.nlme)

sitka.nlme3 <- update(sitka.nlme,
                       corr = corCAR1(0.95, ~Time | tree))
summary(sitka.nlme3)

Fpl <- deriv(~ A + (B-A)/(1 + exp((log(d) - ld50)/th)),
   c("A","B","ld50","th"), function(d, A, B, ld50, th) {})

data(Rabbit)
st <- coef(nls(BPchange ~ Fpl(Dose, A, B, ld50, th),
          start = c(A=25, B=0, ld50=4, th=0.25),
          data = Rabbit))
Rc.nlme <- nlme(BPchange ~ Fpl(Dose, A, B, ld50, th),
    fixed = list(A ~ 1, B ~ 1, ld50 ~ 1, th ~ 1),
    random = A + ld50 ~ 1 | Animal, data = Rabbit,
    subset = Treatment=="Control",
    start = list(fixed=st))
Rm.nlme <- update(Rc.nlme, subset = Treatment=="MDL")

Rc.nlme
Rm.nlme

options(contrasts=c("contr.treatment", "contr.poly"))
c1 <- c(28, 1.6, 4.1, 0.27, 0)
R.nlme1 <- nlme(BPchange ~ Fpl(Dose, A, B, ld50, th),
                fixed = list(A ~ Treatment, B ~ Treatment,
                             ld50 ~ Treatment, th ~ Treatment),
                random =  A + ld50 ~ 1 | Animal/Run, data = Rabbit,
                start = list(fixed=c1[c(1,5,2,5,3,5,4,5)]))
summary(R.nlme1)

R.nlme2 <- update(R.nlme1,
     fixed = list(A ~ 1, B ~ 1, ld50 ~ Treatment, th ~ 1),
     start = list(fixed=c1[c(1:3,5,4)]))
anova(R.nlme2, R.nlme1)
summary(R.nlme2)

if(require(lattice)) {
trellis.device(postscript, file="ch08b.ps", width=8, height=6, pointsize=9)

print(xyplot(BPchange ~ log(Dose) | Animal * Treatment, Rabbit,
   xlab = "log(Dose) of Phenylbiguanide",
   ylab = "Change in blood pressure (mm Hg)",
   subscripts = T, aspect = "xy", panel =
      function(x, y, subscripts) {
         panel.grid()
         panel.xyplot(x, y)
         sp <- spline(x, fitted(R.nlme2)[subscripts])
         panel.xyplot(sp$x, sp$y, type="l")
      }))

} else {

coplot(BPchange ~ log(Dose) | Animal * Treatment, Rabbit,
       show.given=FALSE,
       xlab = "log(Dose) of Phenylbiguanide",
       ylab = "Change in blood pressure (mm Hg)",
       subscripts = TRUE,
       panel = function(x, y, subscripts, ...) {
           points(x, y)
           lines(spline(x, fitted(R.nlme2)[subscripts]))
       })

}

# End of ch08
