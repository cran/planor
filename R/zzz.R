.onAttach <- function(libname, pkgname)
  {
    options(bigmemory.typecast.warning=FALSE)
     options(bigmemory.allow.dimnames=TRUE)
   packageStartupMessage("Loaded planor ", as.character(utils::packageVersion("planor")),"\n")
  }
