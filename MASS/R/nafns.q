napredict <- function(omit, object) UseMethod("napredict")
napredict.omit <- napredict.default <- function(omit, object) object

naresid <- function(omit, object) UseMethod("naresid")
naresid.omit <- naresid.default <- function(omit, object) object
