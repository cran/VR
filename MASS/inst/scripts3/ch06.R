#-*- R -*-

## Script for the Third Edition of `Modern Applied Statistics with S-PLUS'

# Chapter 6   Linear Statistical Models

library(MASS)
postscript(file="ch06.ps", width=8, height=6, pointsize=9)
options(contrasts=c("contr.helmert", "contr.poly"))

# 6.1  A linear regression example

data(whiteside)
if(require(lattice)) { ## can run this is lattice is available

trellis.device(postscript, file="ch06b.ps", width=8, height=6, pointsize=9)
print(xyplot(Gas ~ Temp | Insul, whiteside, panel =
 function(x, y, ...) {
    panel.xyplot(x, y, ...)
    panel.lmline(x, y, ...)
 }, xlab = "Average external temperature (deg. C)",
 ylab = "Gas consumption  (1000 cubic feet)"))
dev.off()

} else {

coplot(Gas ~ Temp | Insul, whiteside, panel =
  function(x, y, ...) {
     points(x, y, ...)
     abline(lm(y ~ x))
  }, xlab = "Average external temperature (deg. C)",
  ylab = "Gas consumption  (1000 cubic feet)")

}

gasB <- lm(Gas ~ Temp, whiteside, subset = Insul=="Before")
gasA <- update(gasB, subset = Insul=="After")

summary(gasB)
summary(gasA)

varB <- deviance(gasB)/gasB$df.resid # direct calculation
varB <- summary(gasB)$sigma^2        # alternative

gasBA <- lm(Gas ~ Insul/Temp - 1, whiteside)
summary(gasBA)

gasQ <- lm(Gas ~ Insul/(Temp + I(Temp^2)) - 1, whiteside)
summary(gasQ)$coef

gasPR <- lm(Gas ~ Insul + Temp, whiteside)
anova(gasPR, gasBA)

oldcon <- options(contrasts = c("contr.treatment", "contr.poly"))
gasBA1 <- lm(Gas ~ Insul*Temp, whiteside)
summary(gasBA1)$coef
options(oldcon)


# 6.2  Model formulae and model matrices

dat <- data.frame(a = factor(rep(1:3, 3)),
                  y = rnorm(9, rep(2:4, 3), 0.1))
obj <- lm(y ~ a, dat)
alf.star <- coef(obj)
alf.star
Ca <- contrasts(dat$a)      # contrast matrix for `a'
alf <- drop(Ca %*% alf.star[-1])
alf
dummy.coef(obj)


N <- factor(Nlevs <- c(0,1,2,4))
contrasts(N)
contrasts(ordered(N))

N2 <- N
contrasts(N2, 2) <- poly(Nlevs, 2)
N2 <- C(N, poly(Nlevs, 2), 2)       # alternative
contrasts(N2)

fractions(ginv(contr.helmert(n = 4)))

Cp <- diag(-1, 4, 5);  Cp[row(Cp) == col(Cp) - 1] <- 1
Cp
fractions(ginv(Cp))


# 6.3  Regression diagnostics

data(hills)
hills.lm <- lm(time ~ dist + climb, hills)
hills.lm
frame()
par(fig=c(0, 0.6, 0, 0.55))
plot(fitted(hills.lm), studres(hills.lm))
abline(h=0, lty=2)
# identify(fitted(hills.lm), studres(hills.lm), row.names(hills))
par(fig=c(0.6, 1, 0, 0.55), pty="s")
qqnorm(studres(hills.lm))
qqline(studres(hills.lm))
par(pty="m")
hills.hat <- lm.influence(hills.lm)$hat
cbind(hills, lev=hills.hat)[hills.hat > 3/35, ]
cbind(hills, pred=predict(hills.lm))["Knock Hill", ]
hills1.lm <- lm(time ~ dist + climb, hills[-18, ])
hills1.lm
lm(time ~ dist + climb, hills[-c(7,18), ])
summary(hills1.lm)
summary(lm(time ~ dist + climb, hills[-18, ], weight=1/dist^2))
lm(time ~ -1 + dist + climb, hills[-18, ], weight=1/dist^2)

#hills <- hills   # make a local copy (needed in S-PLUS >= 5)
hills$ispeed <- hills$time/hills$dist
hills$grad <- hills$climb/hills$dist
hills2.lm <- lm(ispeed ~ grad, hills[-18, ])
hills2.lm
frame()
par(fig=c(0, 0.6, 0, 0.55))
plot(hills$grad[-18], studres(hills2.lm), xlab="grad")
abline(h=0, lty=2)
# identify(hills$grad[-18], studres(hills2.lm), row.names(hills)[-18])
par(fig=c(0.6, 1, 0, 0.55), pty="s")
qqnorm(studres(hills2.lm))
qqline(studres(hills2.lm))
par(pty="m")
hills2.hat <- lm.influence(hills2.lm)$hat
cbind(hills[-18,], lev=hills2.hat)[hills2.hat > 1.8*2/34, ]


# 6.4  Safe prediction

data(wtloss)
quad1 <- lm(Weight ~ Days + I(Days^2), wtloss)
quad2 <- lm(Weight ~ poly(Days, 2), wtloss)
new.x <- data.frame(Days = seq(250, 300, 10),
                    row.names=seq(250, 300, 10))
predict(quad1, newdata=new.x)
predict(quad2, newdata=new.x)
#predict.gam(quad2, newdata=new.x)  # R does not have predict.gam

new.white <- data.frame(Temp = 5*(0:2),
              Insul = factor(rep("After", 3)))
predict(gasBA, new.white) # gives error in earlier S

new.white <- data.frame(Temp = 5*(0:2),
              Insul = factor(rep("After", 3),
                levels=c("After", "Before")))
predict(gasBA, new.white) # wrong in earlier S


# 6.5  Robust and resistant regression

#haveLMR <- (version$major >= 5) | (version$major == 4  && version$minor >=5)
# library(lqs)
data(phones)
phones.lm <- lm(calls ~ year, phones)
attach(phones); plot(year, calls); detach()
abline(phones.lm$coef)
abline(rlm(calls ~ year, phones, maxit=50), lty=2, col=2)
abline(lqs(calls ~ year, phones), lty=3, col=3)
#legend(locator(1), legend=c("least squares", "M-estimate", "LTS", lty=1:3, col=1:3)

summary(lm(calls ~ year, data=phones), cor=F)
summary(rlm(calls ~ year, maxit=50, data=phones), cor=F)
summary(rlm(calls ~ year, scale.est="proposal 2", data=phones), cor=F)
summary(rlm(calls ~ year, data=phones, psi=psi.bisquare), cor=F)

lqs(calls ~ year, data=phones)
lqs(calls ~ year, data=phones, method="lms")
lqs(calls ~ year, data=phones, method="S")

summary(rlm(calls ~ year, data=phones, method="MM"), cor=F)
#if(haveLMR) {
#  phones.lmr <- lmRobMM(calls ~ year, data=phones)
#  print(summary(phones.lmr))
#}

hills.lm
hills1.lm # omitting Knock Hill
rlm(time ~ dist + climb, hills)
summary(rlm(time ~ dist + climb, hills, weights=1/dist^2,
            method="MM"), cor=F)
lqs(time ~ dist + climb, data=hills, nsamp="exact")
summary(hills2.lm) # omitting Knock Hill
summary(rlm(ispeed ~ grad, hills), cor=F)
summary(rlm(ispeed ~ grad, hills, method="MM"), cor=F)
#if(haveLMR)
#  summary(lmRobMM(ispeed ~ grad, data=hills))
lqs(ispeed ~ grad, data=hills)



# 6.6  Bootstrapping linear models

library(boot)
fit <- lm(calls ~ year, data=phones)
ph <- data.frame(phones, res=resid(fit), fitted=fitted(fit))
ph.fun <- function(data, i) {
  d <- data
  d$calls <- d$fitted + d$res[i]
  coef(update(fit, data=d))
}
ph.lm.boot <- boot(ph, ph.fun, R=499)
ph.lm.boot

fit <- rlm(calls ~ year, method="MM", data=phones)
ph <- data.frame(phones, res=resid(fit), fitted=fitted(fit))
# next two lines too slow for testing S-PLUS
#ph.rlm.boot <- boot(ph, ph.fun, R=499)
#ph.rlm.boot


# 6.7  Factorial designs and designed experiments

options(contrasts=c("contr.helmert", "contr.poly"))
data(npk)
npk.aov <- aov(yield ~ block + N*P*K, npk)
npk.aov
summary(npk.aov)
alias(npk.aov)
coef(npk.aov)

options(contrasts=c("contr.treatment", "contr.poly"))
npk.aov1 <- aov(yield ~ block + N + K, npk)
summary.lm(npk.aov1)
se.contrast(npk.aov1, list(N=="0", N=="1"), data=npk)
model.tables(npk.aov1, type="means", se=T)

mp <- c("-","+")
NPK <- expand.grid(N=mp, P=mp, K=mp)
NPK

if(F) {
blocks13 <- fac.design(levels=c(2,2,2),
     factor=list(N=mp, P=mp, K=mp), rep=3, fraction=1/2)
blocks46 <- fac.design(levels=c(2,2,2),
   factor=list(N=mp, P=mp, K=mp), rep=3, fraction=~ -N:P:K)

NPK <- design(block = factor(rep(1:6, rep(4,6))),
   rbind(blocks13, blocks46))
i <- order(runif(6)[NPK$block], runif(24))
NPK <- NPK[i,]  # Randomized

lev <- rep(2,7)
factors <- list(S=mp, D=mp, H=mp, G=mp, R=mp, B=mp, P=mp)
Bike <- fac.design(lev, factors, fraction =
   ~ S:D:G + S:H:R + D:H:B + S:D:H:P)
Bike
replications(~ .^2, data=Bike)
}


# 6.8  An unbalanced four-way layout

data(quine)
attach(quine)
table(Lrn, Age, Sex, Eth)

Means <- tapply(Days, list(Eth, Sex, Age, Lrn), mean)
Vars  <- tapply(Days, list(Eth, Sex, Age, Lrn), var)
SD <- sqrt(Vars)
par(mfrow=c(1,2), pty="s")
plot(Means, Vars, xlab="Cell Means", ylab="Cell Variances")
plot(Means, SD, xlab="Cell Means", ylab="Cell Std Devn.")
detach()

boxcox(Days+1 ~ Eth*Sex*Age*Lrn, data = quine, singular.ok = T,
 lambda = seq(-0.05, 0.45, len = 20))

logtrans(Days ~ Age*Sex*Eth*Lrn, data = quine,
   alpha = seq(0.75, 6.5, len=20), singular.ok = T)

quine.hi <- aov(log(Days + 2.5) ~ .^4, quine)
quine.nxt <- update(quine.hi, . ~ . - Eth:Sex:Age:Lrn)
dropterm(quine.nxt, test="F")
quine.lo <- aov(log(Days+2.5) ~ 1, quine)
addterm(quine.lo, quine.hi, test="F")

quine.stp <- stepAIC(quine.nxt,
   scope = list(upper = ~Eth*Sex*Age*Lrn, lower = ~1),
   trace = F)
quine.stp$anova

dropterm(quine.stp, test="F")
quine.3 <- update(quine.stp, . ~ . - Eth:Age:Lrn)
dropterm(quine.3, test="F")
quine.4 <- update(quine.3, . ~ . - Eth:Age)
dropterm(quine.4, test="F")
quine.5 <- update(quine.4, . ~ . - Age:Lrn)
dropterm(quine.5, test="F")

# 6.9  Predicting computer performance

data(cpus)
boxcox(perf ~ syct+mmin+mmax+cach+chmin+chmax, data=cpus,
       lambda=seq(0, 1, 0.1))

cpus1 <- cpus
attach(cpus)
for(v in names(cpus)[2:7])
  cpus1[[v]] <- cut(cpus[[v]], unique(quantile(cpus[[v]])),
                    include.lowest = T)
detach()
boxcox(perf ~ syct+mmin+mmax+cach+chmin+chmax, data = cpus1,
       lambda = seq(-0.25, 1, 0.1))

set.seed(123)
cpus0 <- cpus1[, 2:8]  # excludes names, authors' predictions
cpus.samp <- sample(1:209, 100)
cpus.lm <- lm(log10(perf) ~ ., data=cpus1[cpus.samp,2:8])
test.cpus <- function(fit)
  sqrt(sum((log10(cpus0[-cpus.samp, "perf"]) -
            predict(fit, cpus0[-cpus.samp,]))^2)/109)
test.cpus(cpus.lm)
cpus.lm2 <- stepAIC(cpus.lm, trace=F)
cpus.lm2$anova
test.cpus(cpus.lm2)

res1 <- log10(cpus1[-cpus.samp, "perf"]) -
              predict(cpus.lm, cpus0[-cpus.samp,])
res2 <- log10(cpus1[-cpus.samp, "perf"]) -
              predict(cpus.lm2, cpus0[-cpus.samp,])
library(ctest)
wilcox.test(res1^2, res2^2, paired=T, alternative="greater")


# 6.10  Multiple comparisons

data(immer)
immer.aov <- aov((Y1+Y2)/2 ~ Var + Loc, data=immer)
summary(immer.aov)

model.tables(immer.aov, type="means", se=T, cterms="Var")

if(F) {
  print(multicomp(immer.aov, plot=T))

  oats1 <- aov(Y ~ N + V + B, data=oats)
  print(summary(oats1))
  print(multicomp(oats1, focus="V"))

  lmat <- matrix(c(0,-1,1,rep(0, 11), 0,0,-1,1, rep(0,10),
                   0,0,0,-1,1,rep(0,9)),,3, dimnames=list(NULL,
                 c("0.2cwt-0.0cwt", "0.4cwt-0.2cwt", "0.6cwt-0.4cwt")))
  print(multicomp(oats1, lmat=lmat, bounds="lower", comparisons="none"))
}


# 6.11  Random and mixed effects

if(F) {
summary(raov(Conc ~ Lab/Bat, data = coop, subset = Spc=="S1"))

#coop <- coop  # make a local copy (needed in S-PLUS >= 5)
is.random(coop) <- T
is.random(coop$Spc) <- F
is.random(coop)
varcomp(Conc ~ Lab/Bat, data=coop, subset = Spc=="S1")
varcomp(Conc ~ Lab/Bat, data=coop, subset = Spc=="S1",
    method = c("winsor", "minque0"))
}

data(oats)
#oats <- oats  # make a local copy: needed in S-PLUS >= 5
oats$Nf <- ordered(oats$N, levels=sort(levels(oats$N)))
oats.aov <- aov(Y ~ Nf*V + Error(B/V), data = oats, qr = T)
summary(oats.aov)
#summary(oats.aov, split=list(Nf=list(L=1, Dev=2:3)))

par(mfrow=c(1,2), pty="s")
plot(fitted(oats.aov[[4]]), studres(oats.aov[[4]]))
abline(h=0, lty=2)
oats.pr <- proj(oats.aov)
qqnorm(oats.pr[[4]][,"Residuals"], ylab="Stratum 4 residuals")
qqline(oats.pr[[4]][,"Residuals"])

par(mfrow=c(1,1), pty="m")
oats.aov2 <- aov(Y ~ N + V + Error(B/V), data = oats, qr = T)
model.tables(oats.aov2, type = "means", se = T)

if(F){
is.random(oats$B) <- T
varcomp(Y ~ N + V + B/V, data = oats)

xyplot(Y ~ EP | No, data = petrol,
    xlab = "ASTM end point (deg. F)",
    ylab = "Yield as a percent of crude",
    panel = function(x, y) {
       m <- sort.list(x)
       panel.grid()
       panel.xyplot(x[m], y[m], type = "b", cex = 0.5)
    })
}

data(petrol)
Petrol <- petrol
names(Petrol)
Petrol[, 2:5] <- scale(as.matrix(Petrol[, 2:5]), scale = F)
pet1.lm <- lm(Y ~ No/EP - 1, Petrol)
matrix(round(coef(pet1.lm),2), 2, 10, byrow = T, dimnames =
   list(c("b0","b1"),levels(Petrol$No)))

pet2.lm <- lm(Y ~ No - 1 + EP, Petrol)
anova(pet2.lm, pet1.lm)

pet3.lm <- lm(Y ~ SG + VP + V10 + EP, Petrol)
anova(pet3.lm, pet2.lm)

library(nlme)
pet3.lme <- lme(Y ~ SG + VP + V10 + EP,
                random = ~ 1 | No, data = Petrol)
summary(pet3.lme)
pet3.lme <- update(pet3.lme, method = "ML")
summary(pet3.lme)
pet4.lme <- update(pet3.lme, fixed = Y ~ V10 + EP)
anova(pet4.lme, pet3.lme)
coef(pet4.lme)
pet5.lme <- update(pet4.lme, random = ~ 1 + EP | No)
anova(pet4.lme, pet5.lme)

data(coop)
lme(Conc ~ 1, random = ~1 | Lab/Bat, data = coop, subset = Spc=="S1")

options(contrasts = c("contr.treatment", "contr.poly"))
oats.lme <- lme(Y ~ N + V, random = ~1 | B/V, data = oats)
summary(oats.lme)

data(Sitka)
sitka.lme <- lme(size ~ treat*ordered(Time),
                 random = ~1 | tree, data=Sitka)
summary(sitka.lme)

Sitka <- Sitka
attach(Sitka)
Sitka$treatslope <- Time * (treat=="ozone")
detach()
sitka.lme2 <- update(sitka.lme,
    fixed = size ~ ordered(Time) + treat + treatslope)
summary(sitka.lme2)
fitted(sitka.lme2, level=0)[1:5]
fitted(sitka.lme2, level=0)[301:305]

sitka.lme <- lme(size ~ treat*ordered(Time), random = ~1 | tree,
                 data = Sitka, corr = corCAR1(, ~Time | tree))
summary(sitka.lme)


# End of ch06
