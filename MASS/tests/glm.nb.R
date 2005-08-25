## glm.nb with weights
library(MASS)
yeast <- data.frame(cbind(numbers = 0:5, fr = c(213, 128, 37, 18, 3, 1)))

attach(yeast)
n <- rep(numbers, fr)

## fitting using glm.nb with weights - wrong results in 7.2-18
yeast2.fit <- glm.nb(numbers~1, link = log, weights=fr)
summary(yeast2.fit)

## fitting extending the vector and using glm.nb - correct result ##
yeast3.fit<-glm.nb(n~1, link = log)
summary(yeast3.fit)

stopifnot(all.equal(deviance(yeast2.fit), deviance(yeast3.fit)))
stopifnot(all.equal(yeast2.fit$theta, yeast3.fit$theta))

detach(yeast)
