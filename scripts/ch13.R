#-*- R -*-

# Chapter 13   Time Series

library(MASS)
postscript(file="ch13.ps", width=8, height=6, pointsize=9)
options(width=65, digits=5)
library(ts)

data(lh)
lh
data(deaths)
deaths
tsp(deaths)
start(deaths)
end(deaths)
frequency(deaths)
#units(deaths)
cycle(deaths)
plot(lh)
data(mdeaths)
data(fdeaths)
plot(ts.union(deaths, mdeaths, fdeaths), lty=c(1,3,4))
aggregate(deaths, 4, sum)
aggregate(deaths, 1, mean)


# 13.1  Second-order summaries

acf(lh)
acf(lh, type="covariance")
acf(deaths)
acf(ts.union(mdeaths, fdeaths))
spectrum(lh)
spectrum(deaths)
par(mfrow=c(2,2))
spectrum(lh)
spectrum(lh, spans=3)
spectrum(lh, spans=c(3,3))
spectrum(lh, spans=c(3,5))

spectrum(deaths)
spectrum(deaths, spans=c(3,3))
spectrum(deaths, spans=c(3,5))
spectrum(deaths, spans=c(5,7))
par(mfrow=c(1,1))
spectrum(deaths)
deaths.spc <- spec.pgram(deaths, taper=0, plot=F)
lines(deaths.spc$freq, deaths.spc$spec, lty=3)

par(mfrow=c(1,2))
cpgram(lh)
cpgram(deaths)
par(mfrow=c(1,1))


# 13.2  ARIMA models

#ts.sim <- arima.sim(list(order=c(1,1,0), ar=0.7), n=200)

acf(lh, type="partial")
acf(deaths, type="partial")
lh.ar1 <- ar(lh, F, 1)
cpgram(lh.ar1$resid, main="AR(1) fit to lh")
lh.ar <- ar(lh, order.max=9)
lh.ar$order
lh.ar$aic
cpgram(lh.ar$resid, main="AR(3) fit to lh")

lh.arima1 <- arima0(lh, order=c(1,0,0), include.mean = T)
arima0.diag(lh.arima1)
lh.arima3 <- arima0(lh, order=c(3,0,0), include.mean = T)
arima0.diag(lh.arima3)
lh.arima11 <- arima0(lh, order=c(1,0,1), include.mean = T)
arima0.diag(lh.arima11)
#lh.fore <- arima.forecast(lh1, n=12, model=lh.arima3$model)
lh.fore <- predict(lh.arima3, n.ahead=12)
#lh.fore$mean <- lh.fore$mean + mean(lh)
#ts.plot(lh, lh.fore$mean, lh.fore$mean+2*lh.fore$std.err,
#        lh.fore$mean-2*lh.fore$std.err)
plot(lh, xlim=c(1, 60))
lines(lh.fore$pred, col="blue")
lines(lh.fore$pred+2*lh.fore$se, lty=2, col="blue")
lines(lh.fore$pred-2*lh.fore$se, lty=2, col="blue")


# 13.3  Seasonality

deaths.stl <- stl(deaths, "periodic")
plot(deaths.stl)
dsd <-  deaths.stl$time.series[, "trend"] + deaths.stl$time.series[, "remainder"]
plot(dsd)
acf(dsd)
acf(dsd, type="partial")
spectrum(dsd, span=c(3,3))
cpgram(dsd)
dsd.ar <- ar(dsd)
dsd.ar$order
dsd.ar$aic
dsd.ar$ar
cpgram(dsd.ar$resid, main="AR(1) residuals")
if(F) {
dsd.rar <- ar.gm(dsd)
dsd.rar$ar
}
deaths.diff <- diff(deaths, 12)
acf(deaths.diff, 30)
acf(deaths.diff, 30, type="partial")
ar(deaths.diff)
# this suggests the seasonal effect is still present.
deaths.arima1 <- arima0(deaths, order=c(2,0,0),
                        list(order=c(0,1,0), period=12))
deaths.arima1
arima0.diag(deaths.arima1, gof.lag=24)
# suggests need a seasonal AR term
deaths.arima2 <- arima0(deaths, order=c(2,0,0),
                        list(order=c(1,0,0), period=12))
deaths.arima2
arima0.diag(deaths.arima2, gof.lag=24)
cpgram(deaths.arima2$resid)
deaths.arima3 <- arima0(deaths, order=c(2,0,0),
                        list(order=c(1,1,0), period=12))
deaths.arima3  # aic not comparable to deaths.arima2
arima0.diag(deaths.arima3, gof.lag=24)
arima0(deaths,order=c(2,0,0), list(order=c(1,0,0), period=12))
arima0(deaths, order=c(2,0,0), list(order=c(2,0,0), period=12))

data(nottem)
nott <- window(nottem, end=c(1936,12))
plot(nott)
nott.stl <- stl(nott, "period")
plot(nott.stl)
nott.stl <- stl(nott, 5)
plot(nott.stl)

par(mfrow=c(1,1))
boxplot(split(nott, cycle(nott)), names.x=month.abb)

nott[110] <- 35
nott.stl <- stl(nott, "period")
nott1 <- nott.stl$time.series[, "trend"] + nott.stl$time.series[, "remainder"]
acf(nott1)
acf(nott1, type="partial")
cpgram(nott1)
ar(nott1)$aic
plot(0:23, ar(nott1)$aic, xlab="order", ylab="AIC", main="AIC for AR(p)")
nott1.ar1 <- arima0(nott1, order=c(1,0,0))
nott1.ar1

nott.fore <- predict(nott1.ar1, n.ahead=36)
nott.fore$pred <- nott.fore$pred +
    as.vector(nott.stl$time.series[1:36, "seasonal"])
plot(window(nottem, 1937), lty=3, ylim=c(30, 70))
lines(nott.fore$pred, col="blue")
lines(nott.fore$pred + 2**nott.fore$se, lty=2, col="blue")
lines(nott.fore$pred - 2**nott.fore$se, lty=2, col="blue")

# nott1.fore <- arima.forecast(nott1, n=36,
#                              model=nott1.ar1$model)
# nott1.fore$mean <- nott1.fore$mean + mean(nott.stl$rem) +
#                         as.vector(nott.stl$sea[1:36])
# ts.plot(window(nottem, 1937), nott1.fore$mean,
#         nott1.fore$mean+2*nott1.fore$std.err,
#         nott1.fore$mean-2*nott1.fore$std.err, lty=c(3,1,2,2))
title("via Seasonal Decomposition")


acf(diff(nott,12), 30)
acf(diff(nott,12), 30, type="partial")
cpgram(diff(nott,12))
nott.arima1 <- arima0(nott,order=c(1,0,0), list(order=c(2,1,0), period=12))
nott.arima1
arima0.diag(nott.arima1, gof.lag=24)
nott.fore <- predict(nott.arima1, n.ahead=36)
plot(window(nottem, 1937), lty=3, ylim=c(30, 70))
lines(nott.fore$pred, col="blue")
lines(nott.fore$pred + 2**nott.fore$se, lty=2, col="blue")
lines(nott.fore$pred - 2**nott.fore$se, lty=2, col="blue")

# nott.fore <- arima.forecast(nott, n=36,
#     model=nott.arima1$model)
# ts.plot(window(nottem, 1937), nott.fore$mean,
#     nott.fore$mean+2*nott.fore$std.err,
#     nott.fore$mean-2*nott.fore$std.err, lty=c(3,1,2,2))
title("via Seasonal ARIMA model")


# 13.6  Regression with autocorrelated errors

data(beav1); data(beav2)
# beav1 <- beav1; beav2 <- beav2
attach(beav1)
beav1$hours <- 24*(day-346) + trunc(time/100) + (time%%100)/60
detach()
attach(beav2)
beav2$hours <- 24*(day-307) + trunc(time/100) + (time%%100)/60
detach()
par(mfrow=c(2,2))
plot(beav1$hours, beav1$temp, type="l", xlab="time",
     ylab="temperature", main="Beaver 1")
usr <- par("usr"); usr[3:4] <- c(-0.2, 8); par(usr=usr)
lines(beav1$hours, beav1$activ, type="s", lty=2)
plot(beav2$hours, beav2$temp, type="l", xlab="time",
     ylab="temperature", main="Beaver 2")
usr <- par("usr"); usr[3:4] <- c(-0.2, 8); par(usr=usr)
lines(beav2$hours, beav2$activ, type="s", lty=2)

attach(beav2)
temp <- ts(temp, start=8+2/3, frequency=6)
activ <- ts(activ, start=8+2/3, frequency=6)
acf(temp[activ==0]); acf(temp[activ==1]) # also look at PACFs
ar(temp[activ==0]); ar(temp[activ==1])
par(mfrow=c(1,1))

arima0(temp, order=c(1,0,0))
arima0(temp, order=c(1,0,0), xreg=activ)

dreg <- cbind(sin=sin(2*pi*hours/24), cos=cos(2*pi*hours/24))
arima0(temp, order=c(1,0,0), xreg=cbind(active=activ,dreg))


alpha <- 0.8255
stemp <- as.vector(temp - alpha*lag(temp, -1))
X <- cbind(1, activ); sX <- X - alpha*lag(X, -1)
beav2.ls <- lm(stemp ~ -1 + sX)
beav2.sls <- summary(beav2.ls)
beav2.sls
sqrt(t(c(1,1)) %*% beav2.sls$cov %*% c(1,1)) * beav2.sls$sigma
plot(hours[-1], residuals(beav2.ls))
detach(); rm(temp, activ)

# lme3
library(lme)
beav2.gls <- gls(temp ~ activ, data=beav2,  corr=corAR1(0.8), method="ML")
summary(beav2.gls)
summary(update(beav2.gls, subset=6:100))

attach(beav1)
temp <- ts(c(temp[1:82], NA, temp[83:114]), start=9.5,
            frequency=6)
activ <- ts(c(activ[1:82], NA, activ[83:114]), start=9.5,
             frequency=6)
acf(temp[1:53])
acf(temp[1:53], type="partial")
ar(temp[1:53])

act <- c(rep(0, 10), activ)
X <- cbind(1, act=act[11:125], act1 = act[10:124],
          act2 = act[9:123], act3 = act[8:122])
#arima0(temp, xreg=X, order=c(1,0,0)) ## currently fails due to NAs

alpha <- 0.80
stemp <- as.vector(temp - alpha*lag(temp, -1))
sX <- X[-1, ] - alpha * X[-115,]
beav1.ls <- lm(stemp ~ -1 + sX, na.action=na.omit)
summary(beav1.ls, cor=F)
detach(); rm(temp, activ)



# End of ch13

