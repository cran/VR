.noGenerics <- TRUE

.onUnload <- function(libpath)
    library.dynam.unload("spatial", libpath)
