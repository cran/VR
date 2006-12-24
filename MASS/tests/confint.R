library(MASS)

PropCI <- function(x, n, conf = 0.95)
{
    DF <- data.frame(y = x / n, weights = n)
    mod <- glm(y ~ 1, weights = weights, family = binomial(), data = DF)
    plogis(confint(mod, level = conf))
}

PropCI(14, 35)
## had scope error prior to 7.2-31
