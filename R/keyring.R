#---------------------------------------------------------------------------
# CLASS "keyring" and its METHODS
#---------------------------------------------------------------------------
# No roxygen documentation maintained, see directly the .Rd file
#---------------------------------------------------------------------------
# An S4 class to represent a collection of key matrices associated with the same prime, with each key matrix a possible solution to the same model and estimate specifications. An object of class \code{\linkS4class{listofkeyrings}} is a list of \code{\linkS4class{keyring}} objects associated with the different primes involved in a given factorial design problem
# @slot main a list of \code{\linkS4class{keymatrix}} objects associated with the same prime \code{p} and with the same factors
# @slot p a prime number
# @slot LIB common row and column names of the key matrices
# @slot pseudo.info a dataframe with one row per pseudofactor and with columns containing information on the factors (or pseudofactors) associated with the columns of the key matrices (e.g. 'parent', 'nlev', 'block' 'ordered', 'model', 'basic'; see the description of the class \code{\linkS4class{designfactors}} 
#---------------------------------------------------------------------------
#  Methods of "keyring": show, summary
#---------------------------------------------------------------------------
setClass("keyring",
         contains=c("list"),
         representation(p="numeric",
                        pseudo.info="data.frame",
                        LIB="list"))

##--------------------------------------------------------------------------
# "show.keyring" help description in roxygen syntax
#' Print an object of class \code{\linkS4class{keyring}}
#'
#' @aliases show,keyring-method
#' @name show.keyring
#' @aliases show.keyring
#' @title  Print  a 'keyring' object
#'   @param object an object of class \code{\linkS4class{keyring}}
#' @return
#' \sQuote{show} returns an invisible \sQuote{NULL}.
#' @author H. Monod, and al.
#' @seealso  \code{\link{pick.listofkeyrings}} or method "["
#' @note
#' - The number of rows and columns of the matrices that are printed
#' are limited by the option \code{planor.max.print}
#'
#' - Non visible slot is: \code{pseudo.info}
#' @examples
#' K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#'    model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#'    nunits=3^3, base=~A+B+C, max.sol=2)
#' ## The method will now be used for automatic printing of a component of K0
#' K0[[1]]
#' show(K0[[1]]) ## idem
#' print(K0[[1]]) ## idem
#' @keywords  design
#' @export
# "End show.keyring" help description in roxygen syntax
# ---------------------------------------------
show.keyring <- function(object){
  cat("An object of class keyring\n")
  cat(" Number of solutions:", length(object),
      "for prime", object@p, "\n\n")
  
  lapply(object, printgmat)
  
  invisible()
} ## fin show.keyring


# --------------------------------------
# "show" method help description in roxygen syntax
# --------------------------------------
#' Method show for "keyring"
#'
#' @title  Method alias for "keyring"
#'
#' @name show-method.keyring
#' @aliases show-method.keyring
setMethod("show", signature(object="keyring"),
          definition=show.keyring)


##--------------------------------------------------------------------------
# "summary.keyring" help description in roxygen syntax
#' Summarizes the design properties from a \code{\linkS4class{keyring}} object, by
#' printing the summary of each of its key matrices
#'
#' @aliases summary,keyring-method
#' @name summary.keyring
#' @aliases summary.keyring
#' @title  Summarize  a 'keyring' object.
#'   @param object an object of class \code{\linkS4class{keyring}}
#'   @param \ldots ignored
#' @note The number of rows and columns of the matrices that are printed
#' are limited by the option \code{planor.max.print}
#' @return
#'    Does not return anything in the present version.
#' @author H. Monod, and al.
#' @seealso  \code{\link{summary.keymatrix}} and the class \code{\linkS4class{keyring}}
#' @examples
#' K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#'    model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#'    nunits=3^3, base=~A+B+C, max.sol=2)
#' summary(K0[[1]])
#' @keywords  design
#' @export
# "End summary.keyring" help description in roxygen syntax
# ---------------------------------------------

summary.keyring <- function(object, show="tbw", save ="kw", ...){
  ## NOTE: the formal argument list "(object, ...)" is
  ## required to be compatible with the generic function
  ## "summary" in R;

  ## Is some display required?
  isshow <-  (length(show) >0 && show != "" &&
    grepl("[d,t,b,w]", show, ignore.case=TRUE)) 
 
  ## Pseudofactors
  pseudo.info <- object@pseudo.info
  FACTtpf <- pseudo.info$parent
  NIVtpf <- pseudo.info$nlev
  LIBtpf <- rownames(pseudo.info)
  BLOCKtpf <- pseudo.info$block
  ## units factors
  PVuqf <- unique(pseudo.info$nlev)
  PVuqf <- PVuqf[order(PVuqf)]
  Nuqf <- length(PVuqf)

  p.k <- object@p
  if (isshow)
    printpower(p.k)
    ## Loop on the key matrices
    nsol <- length(object)
    sortie <- list()   


      for(l in seq_len(nsol)){
        if (isshow)
          cat(paste("--- Solution ", l, " ---\n\n"))
       sortie[[l]] <- summary.keymatrix(object=object[[l]],
                        lib=colnames(object[[l]]),
                        fact=FACTtpf[ NIVtpf == p.k ],
                        block=BLOCKtpf[ NIVtpf == p.k ],
                          show, save)
    } ## fin l
  if ( length(sortie) >0 ) {
    ## some outut are required
      names(sortie) <- paste("Solution", seq_len(nsol))
    return(invisible(sortie))
    }    else {      return(invisible())}

} ## fin summary.keyring

# --------------------------------------
# "summary" method help description in roxygen syntax
# --------------------------------------
#' Method summary for "keyring"
#'
#' @title  Method summary for "keyring"
#'
#' @name summary-method.keyring
#' @aliases summary-method.keyring
setMethod("summary", signature(object="keyring"),
          definition=summary.keyring)

