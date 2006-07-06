.noGenerics <- TRUE

.onAttach <- function(...) require("stats")

.onUnload <- function(libpath)
    library.dynam.unload("MASS", libpath)
