.First.lib <- function(lib, pkg) {
  if(version$major==0 && version$minor < 64)
    stop("This version for R 0.64 or later")
  library.dynam("spatial", pkg, lib)
  assign(".sp.lib.name", file.path(lib, pkg, "data"),
         pos=match("package:spatial", search()))
}

