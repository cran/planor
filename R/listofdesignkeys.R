 
#---------------------------------------------------------------------------
# CLASS "listofdesignkeys" and its METHODS
#---------------------------------------------------------------------------
# No roxygen documentation, see directly the .Rd file
#---------------------------------------------------------------------------
# listofdesignkeys : an S4 class, typically an output from planor.designkey when the research is recursive
#         main: a list of design-key solutions; each component
#                  of main is a whole solution list across the different primes. It is an object of class \code{\linkS4class{designkey}}
#         factors: the 'designfactors' object that defines the factors
#         model: the list of components of type c(model,estimate)
#                containing the model and estimate specifications
#  Methods of "listofdesignkeys": "[" (or pick), planor.design, summary, show, alias
#---------------------------------------------------------------------------
setClass("listofdesignkeys",
         contains=c("list"),
         representation(factors="designfactors",
                        model="list",
                        nunits="numeric"))
#---------------------------------------------------------------------------
# "pick.listofdesignkeys" help description in roxygen syntax
#'   Extract a single \code{\linkS4class{designkey}}
#'   object (with one key matrix per prime)
#'  from an object of class \code{\linkS4class{listofdesignkeys}}
#'
#' @name pick.listofdesignkeys
#' @aliases pick.listofdesignkeys
#' @aliases pick,listofdesignkeys-method
#' @aliases [,listofdesignkeys,ANY,ANY,ANY-method
#' @title Method to extract a design key from a listofdesignkeys object
#'
#' @param keys an object of class \code{\linkS4class{listofdesignkeys}}
#' @param selection an integer equal to the position of the required solution
#' @return  An object of class \code{\linkS4class{designkey}}, which contains  the selected design
#' @note  \code{K <- pick.listofdesignkeys(K0,1)} can be also be written
#' \code{K <- pick(K0,1)} or more simply \code{K <- K0[1]}
#' @keywords design
#' @author H. Monod, and al.
#' @seealso The classes \code{\linkS4class{listofdesignkeys}},  \code{\linkS4class{designkey}}
#' @examples
#' F2 <- planor.factors( factors=c("R","C","U","A","B1","B2"), nlevels=c(3,2,2,3,2,2) )
#' M2 <- planor.model( model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2 )
#' K2 <- planor.designkey(factors=F2, model=M2, nunits=12,
#'                       base=~R+C+U, max.sol=2)
#' K2.1 <- pick(K2,1)
#' K2.1 <- K2[1] ## Another way of extracting ([ is synonym of pick)
#' @export
# End "pick.listofdesignkeys" help description in roxygen syntax
# ------------------------------------------------
pick.listofdesignkeys <- function(keys,selection){
  if(getOption("verbose")){
    cat( "Extraction of a design key from an object of class listofdesignkeys\n" )
  }
  if (selection > length(keys)) {
    stop( paste("The selection argument must be smaller than ", length(keys)))
  }



  pickdesign <- keys[[selection]]
  return(pickdesign)
}

# --------------------------------------
# "pick" method help description in roxygen syntax
# --------------------------------------
#' Method pick for "listofdesignkeys"
#'
#' @title  Method pick for "listofdesignkeys"
#'
#' @name pick-method.listofdesignkeys
#' @aliases pick-method.listofdesignkeys
setMethod("pick", signature(keys="listofdesignkeys"),
          definition=pick.listofdesignkeys)

 ##------------------------------------------------------------------------
##  Method to return the \code{\linkS4class{designkey}} object
##  of  the \code{index} solution
##   from a \code{\linkS4class{listofdesignkeys}} object.
##
##  @name [,listofdesignkeys,ANY,ANY,ANY-method
##  @title Method "[" for listofdesignkeys
# --------------------------------------
setMethod("[", "listofdesignkeys",
          definition=function(x,i,j,...,drop){
            if (missing(j))
              x <- pick.listofdesignkeys(x, c(i,...))
            else
              x <- pick.listofdesignkeys(x, c(i,j,...))
            x
          })


 ##------------------------------------------------------------------------
## "planor.design.listofdesignkeys" help description in roxygen syntax
## ---------------------------------------------------------------
#' Build one design from an object of class
#'  \code{\linkS4class{listofdesignkeys}}
#'
#' @title Build a design from a 'listofdesignkeys' object
#'
#' @aliases planor.design,listofdesignkeys-method
#' @name planor.design.listofdesignkeys
#' @aliases planor.design.listofdesignkeys
#' @param key an object of class \code{\linkS4class{listofdesignkeys}}
#' @param selection integer to select the solution
#' @param randomize an optional formula. When set, the final designs are randomized according to it.
#' @return  An object of class  \code{\linkS4class{planordesign}},
#' which contains the design built from the selected key matrices
#' @note Restricted to giving a single design
#' @keywords design
#' @author H. Monod, and al.
#' @seealso Classes  \code{\linkS4class{planordesign}},
#'  \code{\linkS4class{listofdesignkeys}}
#' @examples
#' K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"), nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2, nunits=12, base=~R+C+U, max.sol=2)
#' P0 <- planor.design(key=K0, select=1)
#' @export
# End  "planor.design.listofdesignkeys" help description in roxygen syntax
# -----------------------------------------------
planor.design.listofdesignkeys <- function(key, randomize=NULL, selection=1, ...){
    selected <- pick.listofdesignkeys(key,selection)
    OUT <- planor.design.designkey(selected, randomize, ...)
    return(OUT)
}
# --------------------------------------
# "planor.design" method help description in roxygen syntax
# --------------------------------------
#' Method planor.design for "listofdesignkeys"
#'
#' @title  Method planor.design for "listofdesignkeys"
#'
#' @name planor.design-method.listofdesignkeys
#' @aliases planor.design-method.listofdesignkeys
setMethod("planor.design", signature(key="listofdesignkeys"),
          definition=planor.design.listofdesignkeys)


##--------------------------------------------------------------------------
# "summary.listofdesignkeys" help description in roxygen syntax
#' Summarizes the design properties of a \code{\linkS4class{listofdesignkeys}} object, by
#' printing the summary of each key matrix in each design key (design key matrix, confounding
#' and aliasing relationships)
#'
#' @aliases summary,listofdesignkeys-method
#' @name summary.listofdesignkeys
#' @aliases summary.listofdesignkeys
#' @title  Summarize  a 'listofdesignkeys' object.
#'   @param object an object of class \code{\linkS4class{listofdesignkeys}}
#'   @param \ldots ignored
#' @return
#'    Does not return anything in the present version.
#' @author H. Monod, and al.
#' @note The number of rows and columns of the matrices that are printed
#' are limited by the option \code{planor.max.print}
#' @seealso  \code{\link{summary.designkey}}, \code{\link{summary.keymatrix}} and the class \code{\linkS4class{listofdesignkeys}}
#' @examples
#' K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"), nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2, nunits=12, base=~R+C+U, max.sol=2)
#' print(summary(K0))
#' @keywords  design
#' @export
# "End summary.listofdesignkeys" help description in roxygen syntax
# ---------------------------------------------

summary.listofdesignkeys <- function(object, show= "tbw", save="kw", ...){
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
  PVuqf <- unique(pseudo.info$nlev)
  PVuqf <- PVuqf[order(PVuqf)]
  Nuqf <- length(PVuqf)

## Loop on the solutions
  nsol <- length(object)
  sortie <- list()
  for(l in seq_len(nsol)){
    if (isshow)
    cat("\n********** Solution", l, "**********\n")
    ## Loop on the distinct prime numbers
    sortie[[l]] <- list()
    for(k in seq_len(Nuqf)){
      p.k <- PVuqf[k]
      if (isshow)
        cat(paste("--- Solution", l, "for prime", p.k, " ---\n\n"))
      retour <- summary.keymatrix(object=object[[l]][[k]],
                        fact=FACTtpf[ NIVtpf == p.k ],
                        block=BLOCKtpf[ NIVtpf == p.k ],
                        show, save)
      if (issave) {
        sortie[[l]][[k]] <- retour
        names(sortie[[l]])[k] <- paste("Solution", l, "for prime", p.k)
      }
    } ## fin k
  } ## fin l
  if ( issave) {
    names(sortie) <- paste("Solution", seq_len(nsol))
    return(invisible(sortie))
  }  else    return(invisible())

} ## fin summary.listofdesignkeys


# --------------------------------------
# "summary" method help description in roxygen syntax
# --------------------------------------
#' Method summary for "listofdesignkeys"
#'
#' @title  Method summary for "listofdesignkeys"
#'
#' @name summary-method.listofdesignkeys
#' @aliases summary-method.listofdesignkeys
setMethod("summary", signature(object="listofdesignkeys"),
          definition=summary.listofdesignkeys)

##--------------------------------------------------------------------------
# "show.listofdesignkeys" help description in roxygen syntax
#' Print the design key matrices of an object of class \code{\linkS4class{listofdesignkeys}}
#'
#' @aliases show,listofdesignkeys-method
#' @name show.listofdesignkeys
#' @aliases show.listofdesignkeys
#' @title  Print the design key matrices of a 'listofdesignkeys' object
#'   @param object an object of class \code{\linkS4class{listofdesignkeys}}
#' @return
#' \sQuote{show} returns an invisible \sQuote{NULL}.
#' @author H. Monod, and al.
#' @seealso The class \code{\linkS4class{listofdesignkeys}}
#' @note The number of rows and columns of the matrices that are printed
#' are limited by the option \code{planor.max.print}
#' @examples
#' K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"), nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2, nunits=12, base=~R+C+U, max.sol=2)
#' K0
#' show(K0) ## idem
#' print(K0) ## idem
#' @keywords  design
#' @export
# "End show.listofdesignkeys" help description in roxygen syntax
# ---------------------------------------------


show.listofdesignkeys <- function(object){
  ## NOTE: the formal argument list "(object)" is
  ## required to be compatible with the generic function
  ##  in R;
  ## units factors
  cat("An object of class listofdesignkeys\n")
  pseudo.info <- object@factors@pseudo.info
  PVuqf <- unique(pseudo.info$nlev)
  PVuqf <- PVuqf[order(PVuqf)]
  Nuqf <- length(PVuqf)

## Loop on the solutions
  nsol <- length(object)
  for(l in seq_len(nsol)){
    cat("\n********** Solution", l, "**********\n")
    ## Loop on the distinct prime numbers
    for(k in seq_len(Nuqf)){
      p.k <- PVuqf[k]
      cat(paste("--- Solution", l, "for prime", p.k, " ---\n\n"))
      printgmat(object[[l]][[k]])
    }
  }

  invisible()
} ## fin show.listofdesignkeys
# --------------------------------------
# "show" method help description in roxygen syntax
# --------------------------------------
#' Method show for "listofdesignkeys"
#'
#' @title  Method alias for "listofdesignkeys"
#'
#' @name show-method.listofdesignkeys
#' @aliases show-method.listofdesignkeys
setMethod("show", signature(object="listofdesignkeys"),
          definition=show.listofdesignkeys)

##--------------------------------------------------------------------------
# "alias.listofdesignkeys" help description in roxygen syntax
#' Summarize the design properties from a \code{\linkS4class{listofdesignkeys}} object.
#'
#' @aliases alias,listofdesignkeys-method
#' @name alias.listofdesignkeys
#' @aliases alias.listofdesignkeys
#' @title  Calculates aliases of a 'listofdesignkeys' object.
#'   @param object an object of class \code{\linkS4class{listofdesignkeys}}
#'   @param model an optional model formula (by default the first model in object)
#'   @param \ldots ignored
#' @return
#'    To see FUNCTION NOT YET IMPLEMENTED
#' @author H. Monod, and al.
#' @seealso  \code{\link{alias.designkey}} and the class \code{\linkS4class{listofdesignkeys}}
#' @keywords  design
#' @export
# "End alias.listofdesignkeys" help description in roxygen syntax
# ---------------------------------------------
alias.listofdesignkeys  <- function(object, model, ...){
  stop("NOT YET IMPLEMENTED\n")
  # VOIR
}
# --------------------------------------
# "alias" method help description in roxygen syntax
# --------------------------------------
#' Method alias for "listofdesignkeys"
#'
#' @title  Method alias for "listofdesignkeys"
#'
#' @name alias-method.listofdesignkeys
#' @aliases alias-method.listofdesignkeys
setMethod("alias", signature(object="listofdesignkeys"),
          definition=alias.listofdesignkeys)
