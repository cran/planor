#---------------------------------------------------------------------------
# CLASS "listofkeyrings" and its METHODS
#---------------------------------------------------------------------------
# No roxygen documentation maintained, see directly the .Rd file
#---------------------------------------------------------------------------
# An S4 class to store design key solutions when there is only one prime involved or when the solutions are independent between primes
# @slot main a list of \code{\linkS4class{keyring}} objects associated with the different primes involved in the design under construction
# @slot factors a \code{\linkS4class{designfactors}} object
# @slot model a list of components of type c(model,estimate) containing the model and estimate specifications
#---------------------------------------------------------------------------
# Methods of "listofkeyrings": "[" (or pick), planor.design, summary, show, alias
#---------------------------------------------------------------------------
setClass("listofkeyrings",
         contains=c("list"),
         representation(factors="designfactors",
                        model="list",
                        nunits="numeric"))
#---------------------------------------------------------------------------

# "pick.listofkeyrings" help description in roxygen syntax
#'   Extract a single \code{\linkS4class{designkey}}
#'   object (with one key matrix per prime)
#'  from an object of class \code{\linkS4class{listofkeyrings}}
#'
#' @name pick.listofkeyrings
#' @aliases pick.listofkeyrings
#' @aliases pick,listofkeyrings-method
#' @aliases [,listofkeyrings,ANY,ANY,ANY-method
#' @title Method to extract a design key from a listofkeyrings object
#'
#' @param keys an object of class \code{\linkS4class{listofkeyrings}}
#' @param selection index vector to select the key matrix for each prime
#' @return  An object of class \code{\linkS4class{designkey}}, which contains  the selected design
#' @note  \code{K <- pick.listofkeyrings(K0,1)} can be also be written
#' \code{K <- pick(K0,1)} or more simply \code{K <- K0[1]}
#' @keywords design
#' @author H. Monod, and al.
#' @seealso The classes \code{\linkS4class{listofkeyrings}},  \code{\linkS4class{designkey}}
#' @examples
#' K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5), model=~block+(A+B+C+D)^2, estimate=~A+B+C+D, nunits=3^3, base=~A+B+C, max.sol=2)
#' K0.1 <- pick.listofkeyrings(K0,1)
#' K0.1 <- pick(K0,1)## Another way of extracting (pick is a method of the class listofkeyrings)
#' K0.1 <- K0[1] ## Another way of extracting ([ is synonym of pick)
#' @export
# End "pick.listofkeyrings" help description in roxygen syntax
# ------------------------------------------------
pick.listofkeyrings <- function(keys,selection){
  if(getOption("verbose")){
    cat( "Extraction of a design key from an object of class listofkeyrings\n")
  }
  Nuqf <- length(keys)
  if ((lgsel <- length(selection)) != Nuqf) {
    stop(
         paste("The length of the selection vector,",
               lgsel,
               ", is different from the number of primes,", Nuqf))
  }

  pickdesign <- vector("list",length=Nuqf)
  names(pickdesign) <- names(keys)

  for(k in seq_len(Nuqf)){
    pickdesign[[k]] <- keys[[k]]@.Data[[selection[k]]]
  }
  pickdesign <- new("designkey", .Data=pickdesign,
                    factors=keys@factors,
                    nunits=keys@nunits,
                    model=keys@model,
                    recursive=FALSE)
  return(pickdesign)
}
# ------------------------------------------------
##  Method to return the \code{\linkS4class{designkey}} object
##  of  index \code{i} for the  of prime, and index \code{j} for the second value, etc, \ldots,
##   from a \code{\linkS4class{listofkeyrings}} object.
##
##  @name [,listofkeyrings,ANY,ANY,ANY-method
##  @title Method "[" for listofkeyrings
# End "[" method help description in roxygen syntax
# --------------------------------------
setMethod("[", "listofkeyrings",
          definition=function(x,i,j,...,drop){
            if (missing(j))
              x <- pick.listofkeyrings(x, c(i,...))
            else
              x <- pick.listofkeyrings(x, c(i,j,...))
            x
          })
# --------------------------------------
# "pick" method help description in roxygen syntax
# --------------------------------------
#' Method pick for "listofkeyrings"
#'
#' @title  Method pick for "listofkeyrings"
#'
#' @name pick-method.listofkeyrings
#' @aliases pick-method.listofkeyrings
setMethod("pick", signature(keys="listofkeyrings"),
          definition=pick.listofkeyrings)

##------------------------------------------------------------------------
## "planor.design.listofkeyrings" help description in roxygen syntax
## ---------------------------------------------------------------
#' Build one design from an object of class \code{\linkS4class{listofkeyrings}}
#'
#' @title Build a design from a 'listofkeyrings' object
#'
#' @aliases planor.design,listofkeyrings-method
#' @name planor.design.listofkeyrings
#' @aliases planor.design.listofkeyrings
#' @param key an object of class \code{\linkS4class{listofkeyrings}}
#' @param selection index vector to select the key matrix for each prime
#' @param randomize an optional formula. When set, the final designs are randomized according to it.
#' @return  An object of class  \code{\linkS4class{planordesign}},
#' which contains the design built from the selected key matrices
#' @note Restricted to giving a single design
#' @keywords design
#' @author H. Monod, and al.
#' @seealso Classes  \code{\linkS4class{planordesign}},
#'  \code{\linkS4class{listofkeyrings}}
#' @examples
#' K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#'              model=~block+(A+B+C+D)^2, estimate=~A+B+C+D, nunits=3^3,
#'             base=~A+B+C, max.sol=2, verbose=TRUE)
#' P0 <- planor.design.listofkeyrings(key=K0, select=2)
#' P0 <- planor.design(key=K0, select=2) ## Another way (planor.design is a method of the class listofkeyrings)
#' P0.R <- planor.design(key=K0, select=2, randomize=~A+B+C+D) ## Randomize the final designs
#' @export
# End  "planor.design.listofkeyrings" help description in roxygen syntax
# -----------------------------------------------
planor.design.listofkeyrings <- function(key, randomize=NULL, selection, ...){
  if(missing(selection)){ selection <- rep(1,length(key)) }
  selected <- pick.listofkeyrings(key,selection)
  OUT <- planor.design.designkey(key=selected, randomize=randomize, ...)
  return(OUT)
}

# --------------------------------------
# "planor.design" method help description in roxygen syntax
# --------------------------------------
#' Method planor.design for "listofkeyrings"
#'
#' @title  Method planor.design for "listofkeyrings"
#'
#' @name planor.design-method.listofkeyrings
#' @aliases planor.design-method.listofkeyrings
setMethod("planor.design", signature(key="listofkeyrings"),
          definition=planor.design.listofkeyrings)
 ##------------------------------------------------------------------------

##--------------------------------------------------------------------------
# "summary.listofkeyrings" help description in roxygen syntax
#' Summarizes the design properties from a \code{\linkS4class{listofkeyrings}} object, by
#' printing the summary of each key matrix in each keyring
#'
#' @aliases summary,listofkeyrings-method
#' @name summary.listofkeyrings
#' @aliases summary.listofkeyrings
#' @title  Summarize  a 'listofkeyrings' object.
#'   @param object an object of class \code{\linkS4class{listofkeyrings}}
#'   @param \ldots ignored
#' @return
#'    Does not return anything in the present version.
#' @author H. Monod, and al.
#' @seealso  \code{\link{summary.designkey}}, \code{\link{summary.keymatrix}} and the class \code{\linkS4class{listofkeyrings}}
#' @note The number of rows and columns of the matrices that are printed
#' are limited by the option \code{planor.max.print}
#' @examples
#' K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#'    model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#'    nunits=3^3, base=~A+B+C, max.sol=2)
#' print(summary(K0))
#' F2 <- planor.factors( factors=c(LETTERS[1:4], "bloc"), nlevels=c(6,6,4,2,6) )
#' M2 <- planor.model( model=~bloc+(A+B+C+D)^2, estimate=~A+B+C+D )
#' K2 <- planor.designkey(factors=F2, model=M2, nunits=144,
#'                        base=~A+B+D, max.sol=2)
#' print(summary(K2))
#' @keywords  design
#' @export
# "End summary.listofkeyrings" help description in roxygen syntax
# ---------------------------------------------

summary.listofkeyrings <- function(object, show="tbw", save="kw", ...){
  ## NOTE: the formal argument list "(object, ...)" is
  ## required to be compatible with the generic function
  ## "summary" in R;

  ## Is some display required?
  isshow <-  (length(show) >0 && show != "" &&
    grepl("[d,t,b,w]", show, ignore.case=TRUE))
  ## Is some output required?
  issave <-  (length(save) >0 && save != "" &&
    grepl("[k,w]", save, ignore.case=TRUE))

  ## Treatment factors
  object@factors <- object@factors[object@factors@fact.info$model]
  fact.info <- object@factors@fact.info
  Ntf <- nrow(fact.info)
  LIBtf <- rownames(fact.info)
  NIVtf <- fact.info$nlev
  BLOCKtf <- fact.info$block
  ## Pseudofactors
  pseudo.info <- object@factors@pseudo.info
  FACTtpf <- pseudo.info$parent
  NIVtpf <- pseudo.info$nlev
  BLOCKtpf <- pseudo.info$block
  ## units factors

  ##PVuqf <- unique(pseudo.info$nlev)
  PVuqf <- unique(factorize(object@nunits))

  PVuqf <- PVuqf[order(PVuqf)]
  Nuqf <- length(PVuqf)
  ## Loop on the distinct prime numbers
    sortie <- list()
  for(k in seq_len(Nuqf)){
    p.k <- PVuqf[k]
    if (isshow)
      printpower(p.k)
    ## Loop on the key matrices
    nsol <- length(object[[k]])
    sortie[[k]] <- list()
     for(l in seq_len(nsol)){
      if (isshow)
        cat(paste("--- Solution ", l, " for prime ", p.k, " ---\n\n"))
      retour <- summary.keymatrix(object=object[[k]][[l]],
                        fact=FACTtpf[ NIVtpf == p.k ],
                        block=BLOCKtpf[ NIVtpf == p.k ],
                        show, save)
      if (issave) {
        sortie[[k]][[l]]  <- retour
        names(sortie[[k]])[l]  <- paste("Solution", l, "for prime", p.k)
      }
    } ## fin l
  } ## fin k
  if ( issave) {
    names(sortie) <- paste("Solution", seq_len(Nuqf))
     return(invisible(sortie))
  }  else    return(invisible())
} ## fin summary.listofkeyrings

# --------------------------------------
# "summary" method help description in roxygen syntax
# --------------------------------------
#' Method summary for "listofkeyrings"
#'
#' @title  Method summary for "listofkeyrings"
#'
#' @name summary-method.listofkeyrings
#' @aliases summary-method.listofkeyrings
setMethod("summary", signature(object="listofkeyrings"),
          definition=summary.listofkeyrings)


##--------------------------------------------------------------------------
# "show.listofkeyrings" help description in roxygen syntax
#' Print the design key matrices of an object of class \code{\linkS4class{listofkeyrings}}
#'
#' @aliases show,listofkeyrings-method
#' @name show.listofkeyrings
#' @aliases show.listofkeyrings
#' @title  Print the design key matrices of a 'listofkeyrings' object
#'   @param object an object of class \code{\linkS4class{listofkeyrings}}
#' @return
#' \sQuote{show} returns an invisible \sQuote{NULL}.
#' @author H. Monod, and al.
#' @seealso The class \code{\linkS4class{listofkeyrings}}
#' @note The number of rows and columns of the matrices that are printed
#' are limited by the option \code{planor.max.print}
#' @examples
#' K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#'   model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#'   nunits=3^3, base=~A+B+C, max.sol=2)
#' ## The method will now be used for automatic printing of a component of K0
#' K0
#' show(K0) ## idem
#' print(K0) ## idem
#' @keywords  design
#' @export
# "End show.listofkeyrings" help description in roxygen syntax
# ---------------------------------------------


show.listofkeyrings <- function(object){
  ## NOTE: the formal argument list "(object)" is
  ## required to be compatible with the generic function
  ##  in R;
  ## units factors
  cat("An object of class listofkeyrings\n")
  pseudo.info <- object@factors@pseudo.info

  ##PVuqf <- unique(pseudo.info$nlev)
  PVuqf <- unique(factorize(object@nunits))

  PVuqf <- PVuqf[order(PVuqf)]
  Nuqf <- length(PVuqf)

  ## A. Design key matrices
  for(k in seq_len(Nuqf)){
    p.k <- PVuqf[k]
    printpower(p.k)
    nsol <- length(object[[k]])
    for(l in seq_len(nsol)){
      cat(paste("--- Solution ", l, " for prime ", p.k, " ---\n\n"))
      printgmat(object[[k]][[l]])
    }
  }

  invisible()
} ## fin show.listofkeyrings
# --------------------------------------
# "show" method help description in roxygen syntax
# --------------------------------------
#' Method show for "listofkeyrings"
#'
#' @title  Method alias for "listofkeyrings"
#'
#' @name show-method.listofkeyrings
#' @aliases show-method.listofkeyrings
setMethod("show", signature(object="listofkeyrings"),
          definition=show.listofkeyrings)


##--------------------------------------------------------------------------
# "alias.listofkeyrings" help description in roxygen syntax
#' Summarize the design properties from a \code{\linkS4class{listofkeyrings}} object.
#' Return the factors, the model and the number of solutions for each prime.
#'
#' @aliases alias,listofkeyrings-method
#' @name alias.listofkeyrings
#' @aliases alias.listofkeyrings
#' @title  Calculates aliases of a 'listofkeyrings' object.
#'   @param object an object of class \code{\linkS4class{listofkeyrings}}
#'   @param model an optional model formula (by default the first model in object)
#'   @param \ldots ignored
#' @return
#'    A list indexed by the primes p of the object. Each element is a 3-column
#'     matrix with one row per solution for prime p. The columns
#'     give (i) the number of unaliased treatment effecs; (ii) the number of
#'     mutually aliased treatment effects;  (iii) the number of treatment effects
#'     aliased with block effects
#' @author H. Monod, and al.
#' @seealso  \code{\link{alias.designkey}} and the class \code{\linkS4class{listofkeyrings}}
#' @examples
#' K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#'    model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#'    nunits=3^3, base=~A+B+C, max.sol=2)
#' alias(K0)
#' F2 <- planor.factors( factors=c(LETTERS[1:4], "bloc"), nlevels=c(6,6,4,2,6) )
#' M2 <- planor.model( model=~bloc+(A+B+C+D)^2, estimate=~A+B+C+D )
#' K2 <- planor.designkey(factors=F2, model=M2, nunits=144,
#'                        base=~A+B+D, max.sol=2)
#' alias(K2)
#' @keywords  design
#' @export
# "End alias.listofkeyrings" help description in roxygen syntax
# ---------------------------------------------

alias.listofkeyrings <- function(object, model, ...){
  if(missing(model)) model <- object@model[[1]][[1]]
  ## NOTE: the formal argument list "(object, ...)" is
  ## required to be compatible with the generic function
  ## "alias" in R;

  ## Treatment factors
  object@factors <- object@factors[object@factors@fact.info$model]
  fact.info <- object@factors@fact.info
  Ntf <- nrow(fact.info)
  LIBtf <- rownames(fact.info)
  NIVtf <- fact.info$nlev
  BLOCKtf <- fact.info$block
  ## nbrs of pseudofactors per factor
  NPStf <- fact.info$npseudo          #FACT.npseudo
  ## Pseudofactors
  pseudo.info <- object@factors@pseudo.info
  FACTtpf <- pseudo.info$parent
  NIVtpf <- pseudo.info$nlev
  BLOCKtpf <- pseudo.info$block
  ## units factors
  PVuqf <- unique(pseudo.info$nlev)
  PVuqf <- PVuqf[order(PVuqf)]
  Nuqf <- length(PVuqf)

  ## Model  (cf. planor.modelterms)
  Mterms <- attributes(terms(model))
  ModelTerms <- Mterms$factors
  if(Mterms$response > 0){ ModelTerms <- ModelTerms[-1,] }
  if(max(ModelTerms)>1){
    stop("Sorry, nesting in the model formulae not possible") }
  ModelTerms <- cbind(MU=0, ModelTerms)
  mLIBtf <- rownames(ModelTerms)
  ## making of coherent rows between model and estimate
  Nmterms <- ncol(ModelTerms)
  ModelFine <- matrix(0, Ntf, Nmterms,
                      dimnames=list(LIBtf, colnames(ModelTerms)))
  ModelFine[mLIBtf,] <- ModelTerms[mLIBtf,]
  ## take off any all-zero column
  ModelFine <- ModelFine[,apply(ModelFine,2,function(x){sum(x)>0})]
  b.modterms <- as.big.matrix(ModelFine,type="short")

  ## Decomposition into pseudofactors
  ## the function "planor.ineligibleset" can do that, assuming that
  ## we are in an "independent search" case
  b.modset <- planor.ineligibleset(object@factors, b.modterms)

  ## output
  alias.stats <- vector(length=Nuqf,mode="list")
  names(alias.stats) <- PVuqf

  ## Loop on the distinct prime numbers
  for(k in seq_len(Nuqf)){
    p <- PVuqf[k]
    printpower(p)
    ## selection of the adequate rows and columns in the model and estimate matrices
    rows.k <- pseudo.info$nlev == p
    model.cols.k <-
      0 < apply(b.modset[rows.k, , drop=FALSE], 2, sum)
    model.k <- b.modset[rows.k, model.cols.k, drop=FALSE]
    ## alias calculations for prime p
    ## Loop on the key matrices
    nsol <- length(object[[k]])
    ## output
    alias.stats.p <- matrix(NA, nrow=nsol, ncol=3)
    colnames(alias.stats.p) <- c("unaliased","trt.aliased","blc.aliased")
    for(l in seq_len(nsol)){
      cat(paste("--- Solution ", l, " for prime ", p, " ---\n\n"))
      alias.stats.p[l,] <- alias.keymatrix(object=object[[k]][[l]],
                                         model=model.k,
                                         fact=FACTtpf[rows.k],
                                         block=BLOCKtpf[rows.k])
    } ## fin
    cat("--- Synthesis on the aliased treatment effects for prime ",
        p, " ---\n\n")
    print(alias.stats.p)
    alias.stats[[k]] <- alias.stats.p
  } ## fin k

  return(invisible(alias.stats))
} ## fin alias.listofkeyrings
#' Method alias for "listofkeyrings"
#'
#' @title  Method alias for "listofkeyrings"
#'
#' @name alias-method.listofkeyrings
#' @aliases alias-method.listofkeyrings
setMethod("alias", signature(object="listofkeyrings"),
          definition=alias.listofkeyrings)

