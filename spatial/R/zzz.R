.First.lib <- function(lib, pkg)
{
    if(version$major==0 && version$minor < 64)
        stop("This version for R 0.64 or later")
    library.dynam("spatial", pkg, lib)
    assign(".sp.lib.name", paste(lib, pkg, "data", sep="/"),
           pos=match("package:spatial", search()))
}

