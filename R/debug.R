# -----------------------------------------------------

# Permet de voir les objets volumineux
# -----------------------------------------------------
all.object.sizes <- function() {
  res <- NULL
  nomres <- NULL
  for (x in ls(envir=parent.frame())) {
     z<- mget(x, envir=parent.frame(), ifnotfound = list(NULL))
    if (!is.null(zu<-unlist(z, use.names = FALSE))) {
#     res <- c(res, object.size(mget(x, envir=parent.frame(), ifnotfound = list(1))))
     res <- c(res, object.size(zu))
     nomres <- c(nomres, x)
     } else cat(x, "pas trouve\n")
   } # fin for
  names(res) <-nomres
  return(sort(res))
}
