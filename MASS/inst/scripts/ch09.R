#-*- R -*-

# Chapter 9   Modern Regression

library(MASS)
postscript(file="ch09.ps", width=8, height=6, pointsize=9)
set.seed(123)
cpus.samp <- sample(1:209, 100)
data(cpus)


# 9.1  Additive models and scatterplot smoothers}

data(mcycle)
attach(mcycle)
par(mfrow = c(3,2))
plot(times, accel, main="Polynomial regression")
lines(times, fitted(lm(accel ~ poly(times, 3))))
lines(times, fitted(lm(accel ~ poly(times, 6))), lty=3)
legend(40, -100, c("degree=3", "degree=6"), lty=c(1,3),bty="n")
library(splines)
plot(times, accel, main="Natural splines")

lines(times, fitted(lm(accel ~ ns(times, df=5))))
lines(times, fitted(lm(accel ~ ns(times, df=10))), lty=3)
lines(times, fitted(lm(accel ~ ns(times, df=20))), lty=4)
legend(40, -100, c("df=5", "df=10", "df=20"), lty=c(1,3,4),
   bty="n")

library(modreg)
plot(times, accel, main="Smoothing splines")
lines(smooth.spline(times, accel))
plot(times, accel, main="Lowess")
lines(lowess(times, accel))
lines(lowess(times, accel, 0.2), lty=3)
legend(40, -100, c("default span", "f = 0.2"), lty=c(1,3),
    bty="n")

plot(times, accel, main ="ksmooth")
lines(ksmooth(times, accel,"normal", bandwidth=5))
lines(ksmooth(times, accel,"normal", bandwidth=2), lty=3)
legend(40, -100, c("bandwidth=5", "bandwidth=2"), lty=c(1,3),
   bty="n")

plot(times, accel, main ="supsmu")
lines(supsmu(times, accel))
lines(supsmu(times, accel, bass=3), lty=3)
legend(40, -100, c("default", "bass=3"), lty=c(1,3), bty="n")
detach()

data(rock)
#attach(rock)
rock.lm <- lm(log(perm) ~ area + peri + shape, data=rock)
summary(rock.lm)

library(mgcv)
rock.gam <- gam(log(perm) ~ s(area) + s(peri) + s(shape), data=rock)
rock.gam
#summary(rock.gam)
#anova(rock.lm, rock.gam)
par(mfrow=c(2,3), pty="s")
plot(rock.gam, se=TRUE, pages=0)
rock.gam1 <- gam(log(perm) ~ area + peri + s(shape), data=rock)
par(mfrow=c(1,1))
plot(rock.gam1, se=TRUE)
#anova(rock.lm, rock.gam1, rock.gam)
#detach()

cpus0 <- cpus[, 2:8]
for(i in 1:3) cpus0[,i] <- log10(cpus0[,i])
test.cpus <- function(fit)
  sqrt(sum((log10(cpus0[-cpus.samp, "perf"]) -
            predict(fit, cpus0[-cpus.samp,]))^2)/109)
if(F) {
cpus.gam <- gam(log10(perf) ~ ., data=cpus0[cpus.samp, ])
cpus.gam2 <- step.gam(cpus.gam, scope=list(
  "syct"  = ~ 1 + syct + s(syct, 2) + s(syct),
  "mmin"  = ~ 1 + mmin + s(mmin, 2) + s(mmin),
  "mmax"  = ~ 1 + mmax + s(mmax, 2) + s(mmax),
  "cach"  = ~ 1 + cach + s(cach, 2) + s(cach),
  "chmin" = ~ 1 + chmin + s(chmin, 2) + s(chmin),
  "chmax" = ~ 1 + chmax + s(chmax, 2) + s(chmax)
))
print(cpus.gam2$anova, digits=3)
test.cpus(cpus.gam2)
}
cpus.gam <- gam(log10(perf) ~ s(syct) + s(mmin) + s(mmax) +
                s(cach) + s(chmin) + s(chmax), data=cpus0[cpus.samp, ])
cpus.gam
test.cpus(cpus.gam)


if(!exists("bwt")) {
  data(birthwt)
  attach(birthwt)
  race <- factor(race, labels=c("white", "black", "other"))
  ptd <- factor(ptl > 0)
  ftv <- factor(ftv); levels(ftv)[-(1:2)] <- "2+"
  bwt <- data.frame(low=factor(low), age, lwt, race,
	   smoke=(smoke>0), ptd, ht=(ht>0), ui=(ui>0), ftv)
  detach(); rm(race, ptd, ftv)
}

attach(bwt)
age1 <- age*(ftv=="1"); age2 <- age*(ftv=="2+")
birthwt.gam <- gam(low ~ s(age) + s(lwt) + smoke + ptd +
    ht + ui + ftv + s(age1) + s(age2) + smoke:ui, binomial, bwt)
birthwt.gam
table(low, predict(birthwt.gam, bwt) > 0)
if(interactive())plot(birthwt.gam, se=T)
detach()



# 9.2  Projection-pursuit regression

library(modreg)
attach(rock)
rock1 <- data.frame(area=area/10000, peri=peri/10000,
		    shape=shape, perm=perm)
detach()
rock.ppr <- ppr(log(perm) ~ area + peri + shape, data=rock1,
                nterms=2, max.terms=5)
rock.ppr
summary(rock.ppr)

par(mfrow=c(3,2))
plot(rock.ppr)
plot(update(rock.ppr, bass=5))
plot(update(rock.ppr, sm.method="gcv", gcvpen=2))
par(mfrow=c(1,1))

rock.ppr2 <- update(rock.ppr, sm.method="gcv", gcvpen=2)
summary(rock.ppr2)

summary(rock1) # to find the ranges of the variables
Xp <- expand.grid(area=seq(0.1,1.2,0.05),
                  peri=seq(0,0.5,0.02), shape=0.2)
rock.grid <- cbind(Xp,fit=predict(rock.ppr2, Xp))
#wireframe(fit ~ area+peri, rock.grid, screen=list(z=160,x=-60),
#          aspect=c(1,0.5), drape=T)
persp(seq(0.1,1.2,0.05), seq(0,0.5,0.02), matrix(rock.grid$fit,23),
      d=5, theta=-160, phi=30, zlim=c(-1, 15))

# From Chapter 6, for comparisons
cpus1 <- cpus
attach(cpus)
for(v in names(cpus)[2:7])
  cpus1[[v]] <- cut(cpus[[v]], unique(quantile(cpus[[v]])),
                    include.lowest = T)
detach()
cpus.lm <- lm(log10(perf) ~ ., data=cpus1[cpus.samp,2:8])
cpus.lm2 <- stepAIC(cpus.lm, trace=F)
res2 <- log10(cpus1[-cpus.samp, "perf"]) -
              predict(cpus.lm2, cpus1[-cpus.samp,])

cpus.ppr <- ppr(log10(perf) ~ ., data=cpus0[cpus.samp,],
                nterms=2, max.terms=10, bass=5)
cpus.ppr

cpus.ppr <- ppr(log10(perf) ~ ., data=cpus0[cpus.samp,],
               nterms=8, max.terms=10, bass=5)
test.cpus(cpus.ppr)
ppr(log10(perf) ~ ., data=cpus0[cpus.samp,],
    nterms=2, max.terms=10, sm.method="spline")
cpus.ppr2 <- ppr(log10(perf) ~ ., data=cpus0[cpus.samp,],
    nterms=7, max.terms=10, sm.method="spline")
test.cpus(cpus.ppr2)
res3 <- log10(cpus0[-cpus.samp, "perf"]) -
              predict(cpus.ppr, cpus0[-cpus.samp,])
library(ctest)
wilcox.test(res2^2, res3^2, paired=T, alternative="greater")


# 9.3  Response transformation models

library(acepack)
attach(cpus0)
cpus.avas <- avas(cpus0[, 1:6], perf)
plot(log10(perf), cpus.avas$ty)
par(mfrow=c(2,3))
for(i in 1:6) {
  o <- order(cpus0[, i])
  plot(cpus0[o, i], cpus.avas$tx[o, i], type="l",
       xlab=names(cpus0[i]), ylab="")
}
detach()




# 9.4  Neural networks

library(nnet)
attach(rock)
area1 <- area/10000; peri1 <- peri/10000
rock1 <- data.frame(perm, area=area1, peri=peri1, shape)
rock.nn <- nnet(log(perm) ~ area + peri + shape, rock1,
     size=3, decay=1e-3, linout=T, skip=T, maxit=1000, Hess=T)
sum((log(perm) - predict(rock.nn))^2)
detach()
eigen(rock.nn$Hess, T)$values

Xp <- expand.grid(area=seq(0.1,1.2,0.05),
                  peri=seq(0,0.5,0.02), shape=0.2)
rock.grid <- cbind(Xp,fit=predict(rock.nn, Xp))
#wireframe(fit ~ area + peri, rock.grid, screen=list(z=160,x=-60),
#          aspect=c(1,0.5), drape=T)
persp(seq(0.1,1.2,0.05), seq(0,0.5,0.02), matrix(rock.grid$fit,23),
      d=5, theta=-160, phi=30, zlim=c(-1, 15))

attach(cpus0)
cpus1 <-
  data.frame(syct=syct-2, mmin=mmin-3, mmax=mmax-4, cach=cach/256,
             chmin=chmin/100, chmax=chmax/100, perf=perf)
detach()

test.cpus <- function(fit)
  sqrt(sum((log10(cpus1[-cpus.samp, "perf"]) -
           predict(fit, cpus1[-cpus.samp,]))^2)/109)
cpus.nn1 <- nnet(log10(perf) ~ ., cpus1[cpus.samp,], linout=T,
                 skip=T, size=0)
test.cpus(cpus.nn1)

cpus.nn2 <- nnet(log10(perf) ~ ., cpus1[cpus.samp,], linout=T,
                 skip=T, size=4, decay=0.01, maxit=1000)
test.cpus(cpus.nn2)

cpus.nn3 <- nnet(log10(perf) ~ ., cpus1[cpus.samp,], linout=T,
                 skip=T, size=10, decay=0.01, maxit=1000)
test.cpus(cpus.nn3)

cpus.nn4 <- nnet(log10(perf) ~ ., cpus1[cpus.samp,], linout=T,
                 skip=T, size=25, decay=0.01, maxit=1000)
test.cpus(cpus.nn4)


CVnn.cpus <- function(formula, data=cpus1[cpus.samp, ],
    size = c(0, 4, 4, 10, 10),
    lambda = c(0, rep(c(0.003, 0.01), 2)),
    nreps = 5, nifold = 10, ...)
{
  CVnn1 <- function(formula, data, nreps=1, ri,  ...)
  {
    truth <- log10(data$perf)
    res <- numeric(length(truth))
    cat("  fold")
    for (i in sort(unique(ri))) {
      cat(" ", i,  sep="")
      for(rep in 1:nreps) {
        learn <- nnet(formula, data[ri !=i,], trace=F, ...)
        res[ri == i] <- res[ri == i] +
                        predict(learn, data[ri == i,])
      }
    }
    cat("\n")
    sum((truth - res/nreps)^2)
  }
  choice <- numeric(length(lambda))
  ri <- sample(nifold, nrow(data), replace=T)
  for(j in seq(along=lambda)) {
    cat("  size =", size[j], "decay =", lambda[j], "\n")
    choice[j] <- CVnn1(formula, data, nreps=nreps, ri=ri,
                       size=size[j], decay=lambda[j], ...)
    }
  cbind(size=size, decay=lambda, fit=sqrt(choice/100))
}
CVnn.cpus(log10(perf) ~ ., data=cpus1[cpus.samp,],
          linout=T, skip=T, maxit=1000)


# End of ch09
