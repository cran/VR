wt.huber <- function(u, c = 1.345)  ifelse(abs(u) < c, 1, c / abs(u))

.First.lib <- function(lib, pkg) {
  library.dynam("MASS", pkg, lib)
}
