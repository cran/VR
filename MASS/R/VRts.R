diagnostics <- function(x, ...) UseMethod("diagnostics")

VRarima <- function(x, order,
                    seasonal = list(order = c(0, 0, 0), period = NA),
                    n.cond, ...)
{
    Call <- match.call()
    xn <- deparse(substitute(x))
    res <- arima0(x, order, seasonal, ...)
    res$call <- Call
    res$series <- xn
    class(res) <- c("VRarima", "arima0")
    res
}

print.VRarima <-
    function (x, digits = max(3, options("digits")$digits - 3), se = TRUE,
              ...)
{
#    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("Coefficients:\n")
    coef <- round(x$coef, digits = digits)
    if (se && nrow(x$var.coef)) {
        ses <- round(sqrt(diag(x$var.coef)), digits = digits)
        coef <- matrix(coef, 1, dimnames = list(NULL, names(coef)))
        coef <- rbind(coef, "s.e."=ses)
    }
    print.default(coef)
    cat("\nsigma^2 estimated as ", format(x$sigma2, digits = digits),
        ":  log likelihood = ", format(round(x$loglik, 2)), ",  aic = ",
        format(round(x$aic, 2)), "\n", sep = "")
    invisible(x)
}

predict.VRarima <- function(object, n.ahead = 1, se.fit = TRUE, ...)
    predict.arima0(object, n.ahead = n.ahead, se.fit = se.fit, ...)

diagnostics.VRarima <- function(x, ...) arima0.diag(x, ...)
