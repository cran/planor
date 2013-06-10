 
#---------------------------------------------------------------------------
# CLASS designkey and its METHODS
#---------------------------------------------------------------------------
# No roxygen documentation, see directly the .Rd file
#---------------------------------------------------------------------------
# An S4 class to represent a design key solution in the planor package
# main a single design-key solution = a list with one key matrix for each prime. Each one is an objet of class keymatrix
#         factors: the 'designfactors' object that defines the factors
#         model: the list of components of type c(model,estimate)
#                containing the model and estimate specifications
#         recursive: logical, TRUE if the design has been constructed recursively
# Methods of "designkey" : planor.design, summary, show, alias
#---------------------------------------------------------------------------

setClass("designkey",
         contains=c("list"),
         representation(factors="designfactors",
                        model="list",
                        nunits="numeric",
                        recursive="logical"))
##---------------------------------------------------------------------------
## "planor.design.designkey" help description in roxygen syntax
## ----------------------------------------------------
#'  Build the design from a design key matrix
#'
#' @title Build a design from a 'designkey' object
#'
#' @aliases planor.design,designkey-method
#' @name planor.design.designkey
#' @aliases planor.design.designkey
#' @param key an object of class  \code{\linkS4class{designkey}}
#' @param randomize an optional formula. When set, the final design is randomized according to it.
#' @return An object of class \code{\linkS4class{planordesign}},
#' which  contains the design build from the input.
#' @keywords design
#' @author H. Monod, and al.
#' @seealso The classes  \code{\linkS4class{planordesign}}, \code{\linkS4class{designkey} }
#' @examples
#' K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#'   model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#'   nunits=3^3, base=~A+B+C, max.sol=2, verbose=TRUE)
#' P0 <- planor.design.designkey(K0[1])
#' P0 <- planor.design(K0[1]) ## Another way (planor.design is a method of the class designkey)
#' P0.R <- planor.design(K0[1], randomize=~A+B+C+D) ## randomize the final design
#' @export
# End "planor.design.designkey" help description in roxygen syntax
# -------------------------------------------------------
planor.design.designkey <- function(key, randomize=NULL, ...){
    key@factors <- key@factors[key@factors@fact.info$model]
    fact.info <- key@factors@fact.info
    Ntf <- nrow(fact.info)
    LIBtf <- rownames(fact.info)
    NIVtf <- fact.info$nlev
    ## Pseudofactors
    pseudo.info <- key@factors@pseudo.info
    FACTtpf <- pseudo.info$parent
    NIVtpf <- pseudo.info$nlev
    LIBtpf <- rownames(pseudo.info)
    Ntpf <- nrow(pseudo.info)

    ## A. Construction of each sub-design associated with each Sylow subgroup

    ##PVuqf <- unique(pseudo.info$nlev)
    PVuqf <- unique(factorize(key@nunits))

    PVuqf <- PVuqf[order(PVuqf)]
    Nuqf <- length(PVuqf)
    b.pseudodesign.k <- vector("list", length=Nuqf)
    for(k in seq_len(Nuqf)){
        p.k <- PVuqf[k]
        r.k <- nrow(key[[k]])
      b.aux <- crossingBig(rep(p.k,r.k),start=0)
      b.pseudodesign.k[[k]] <- multBigmod(b.aux, key[[k]], p.k)

    }

    ## B. Crossing of the subdesigns
    b.fullpseudodesign <- cross.designs(b.pseudodesign.k)

    ## C. Reordering of the columns by treatment factor
    pseudosorder <- order(NIVtpf, FACTtpf, seq_along(NIVtpf))
    b.fullpseudodesign <- exchangeColBig(b.fullpseudodesign,
                                         pseudosorder,
                                         seq.int(ncol(b.fullpseudodesign)))

    ## D. Calculation of the design with the original treatment factors
    b.back <- big.matrix( Ntpf, Ntf, init=0, type="short")

    for(i in seq_len(Ntpf)){
        select <- (FACTtpf == FACTtpf[i]) & (seq_len(Ntpf) > i)
        b.back[i,FACTtpf[i]] <- prod( NIVtpf[select] )
    }
    b.finaldesign <- multBig(b.fullpseudodesign, b.back)
    ## Transformation de b.finaldesign en non big matrix
    ## car il est factor
    b.finaldesign <- as.data.frame(b.finaldesign[,])

    names(b.finaldesign) <- LIBtf ## noms des colonnes
    for(i in seq_len(Ntf)){
        zz <- factor(b.finaldesign[,i])
        levels(zz) <- (key@factors@levels)[[i]]
        b.finaldesign[i] <- zz ## transfo de la colonne i en facteur
    }

    ## E. Randomization, if required
    if (!is.null(randomize)) {
      b.finaldesign <- planor.randomize(randomize, b.finaldesign, ...)
      if( "InitialUNITS" %in% names(b.finaldesign) ){
        key@factors <- bind(planor.factors( factors="InitialUNITS",
                                           nlevels=nrow(b.finaldesign),
                                           dummy=FALSE),
                            key@factors)
      }
    }

    OUT <- new("planordesign",
               design=b.finaldesign,
               factors=key@factors,
               model=key@model,
               designkey=key@.Data,
               nunits= key@nunits,
               recursive=key@recursive )
    

    return(OUT)
} ## fin planor.design.designkey
##------------------------------------------------------------------------
#' Method planor.design for "designkey"
#'
#' @title  Method planor.design for "designkey"
#'
#' @name planor.design-method.designkey
#' @aliases planor.design-method.designkey
setMethod("planor.design", signature(key="designkey"),
          definition=planor.design.designkey)
##--------------------------------------------------------------------------


# "summary.designkey" help description in roxygen syntax
#' Summarises the design properties of a \code{\linkS4class{designkey}} object, by
#' printing the summary of each of its key matrices (design key matrix, confounding
#' and aliasing relationships)
#'
#' @aliases summary,designkey-method
#' @name summary.designkey
#' @aliases summary.designkey
#' @title  Summarise the design properties of a 'designkey' object
#'   @param object an object of class \code{\linkS4class{designkey}}
#'   @param \ldots ignored
#' @return
#'     A list of matrices, one per prime. Their columns are kernel generators of the key matrices
#' @author H. Monod, and al.
#' @seealso  The class \code{\linkS4class{designkey}}; \code{\link{summary.keymatrix}}
#' @examples
#' K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#'    model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#'    nunits=3^3, base=~A+B+C, max.sol=2)
#' resum <- summary(K0[1])
#' @keywords  design
#' @export
# "End summary.designkey" help description in roxygen syntax
# ---------------------------------------------

summary.designkey <- function(object,show="dtbw", save="k", ...){
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

  ## Storage preparation of the returned results
  Hgen <- vector(mode="list",length=Nuqf)
  ## Loop on the distinct prime numbers
  for(k in seq_len(Nuqf)){
    p.k <- PVuqf[k]
    if (isshow)
      printpower(p.k)
    retour <- summary.keymatrix(object=object[[k]],
                                  fact=FACTtpf[ NIVtpf == p.k ],
                                  block=BLOCKtpf[ NIVtpf == p.k ],
                                   show, save)

    if (issave) {
      Hgen[[k]] <- retour
    }
  } ## fin k

    if (issave) {
      names(Hgen)  <- paste("Prime",  PVuqf[1:Nuqf])
      return(invisible(Hgen))
    } else  return(invisible())
} ## fin summary.designkey

##--------------------------------------------------------------------------
#' Method summary for "designkey"
#'
#' @title  Method summary for "designkey"
#'
#' @name summary-method.designkey
#' @aliases summary-method.designkey
setMethod("summary", signature(object="designkey"),
          definition=summary.designkey)
##--------------------------------------------------------------------------
# "show.designkey" help description in roxygen syntax
#' Print the design key matrices of an object of class \code{\linkS4class{designkey}}
#'
#' @aliases show,designkey-method
#' @name show.designkey
#' @aliases show.designkey
#' @title  Print the design key matrices of a 'designkey' object
#'   @param object an object of class \code{\linkS4class{designkey}}
#' @return
#' \sQuote{show} returns an invisible \sQuote{NULL}.
#' @author H. Monod, and al.
#' @seealso  \code{\link{summary.designkey}} and the class \code{\linkS4class{designkey}}
#' @note
#' - The number of rows and columns of the matrices that are printed
#' are limited by the option \code{planor.max.print}
#'
#' -  Objects of class \code{\linkS4class{designkey}} are displayed automatically is if by a call to 'show'.
#' @examples
#' K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#'    model=~block+(A+B+C+D)^2, estimate=~A+B+C+D ,
#'    nunits=3^3, base=~A+B+C, max.sol=2)
#' ## The method will now be used for automatic printing of a component of K0
#' K0[1]
#' show(K0[1]) ## idem
#' print(K0[1]) ## idem
#' @keywords  design
#' @export
# "End show.designkey" help description in roxygen syntax
# ---------------------------------------------
show.designkey <- function(object){
  ## NOTE: the formal argument list "(object)" is
  ## required to be compatible with the generic function
  ##  in R;
  ## units factors
  cat("An object of class designkey\n")
  keys <- unclass(object)
  primes <- as.integer( names(keys) )
  ##cat("DESIGN KEY MATRICES\n")

  ## A. Design key matrices
  for(k in seq(length(primes))) {
    printpower(primes[k])
    printgmat(keys[[k]])
  }
  invisible()
} ## fin show.designkey

##--------------------------------------------------------------------------
#' Method show for "designkey"
#'
#' @title  Method show for "designkey"
#'
#' @name show-method.designkey
#' @aliases show-method.designkey
setMethod("show", signature(object="designkey"),
          definition=show.designkey)

##--------------------------------------------------------------------------


# "alias.designkey" help description in roxygen syntax
#' Summarise the design properties from a design key matrix.
#' Display the design keys matrices and the factorial effects confounded with the mean.
#'
#' @aliases alias,designkey-method
#' @name alias.designkey
#' @aliases alias.designkey
#' @title Prints the aliases of a 'designkey' object
#'   @param object an object of class \code{\linkS4class{designkey}}
#'   @param model an optional model formula (by default the first model in object)
#'   @param \ldots ignored
#' @return
#'     invisible
#' @author H. Monod, and al.
#' @seealso  The class \code{\linkS4class{designkey}}
#' @examples
#' K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
#'    model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
#'    nunits=3^3, base=~A+B+C, max.sol=2)
#' alias(K0[1])
#' @keywords  design
#' @export
# "End alias.designkey" help description in roxygen syntax
# ---------------------------------------------

alias.designkey <- function(object, model, ...){
  if(missing(model)) model <- object@model[[1]][[1]]
  ## NOTE: the formal argument list "(object, ...)" is
  ## required to be compatible with the generic function
  ## "alias" in R;

  ## Treatment factors
  object@factors <- object@factors[object@factors@fact.info$model]
  fact.info <- object@factors@fact.info
  Ntf <- nrow(fact.info)              #FACT.N
  LIBtf <- rownames(fact.info)
  NIVtf <- fact.info$nlev             #FACT.nlev
  BLOCKtf <- fact.info$block
  ## nbrs of pseudofactors per factor
  NPStf <- fact.info$npseudo          #FACT.npseudo
  ## Pseudofactors
  pseudo.info <- object@factors@pseudo.info
  FACTtpf <- pseudo.info$parent       #PSEUDO.parent
  NIVtpf <- pseudo.info$nlev          #PSEUDO.nlev
  BLOCKtpf <- pseudo.info$block
  ## units factors
  PVuqf <- unique(pseudo.info$nlev)
  PVuqf <- PVuqf[order(PVuqf)]        #PRIMES
  Nuqf <- length(PVuqf)               #UQUASI.N

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
  ## Loop on the distinct primes
  for(k in seq(Nuqf)){
    p <- PVuqf[k]
    printpower(p)
    ## selection of the adequate rows and columns in the model and estimate matrices
    rows.k <- pseudo.info$nlev == p
    modset <- matrix(b.modset[rows.k, ], ncol=ncol(b.modset))
    model.cols.k <-
      0 < apply(modset, 2, sum)
# BUG fixed 30/4/2013    model.k <- modset[rows.k, model.cols.k, drop=FALSE]
    model.k <- modset[, model.cols.k, drop=FALSE]
    ## alias calculations for prime p
    alias.keymatrix(object=object[[k]], model=model.k,
                    fact=FACTtpf[rows.k],
                    block=BLOCKtpf[rows.k])
  } ## fin k

  return(invisible())
} ## fin alias.designkey

##--------------------------------------------------------------------------
#' Method alias for "designkey"
#'
#' @title  Method alias for "designkey"
#'
#' @name alias-method.designkey
#' @aliases alias-method.designkey
setMethod("alias", signature(object="designkey"),
          definition=alias.designkey)
##--------------------------------------------------------------------------
