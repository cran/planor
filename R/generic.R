
##---------------------------------------------------------------------------
# GENERIC FUNCTIONS
##---------------------------------------------------------------------------

# ---------------------------------------------------
# Generic "bind" help description in roxygen syntax
# ---------------------------------------------------
#'  Generic function to bind two objects.
#'
#' @title Bind two objects (generic function)
#'
#' @name bind
#' @aliases bind
#' @param x first object 
#' @param y second object
#' @param ... other arguments
#' @return An object where the input are binded.
#' @keywords design
#' @author H. Monod, and al.
#' @seealso Method \code{\link{bind.designfactors}}
#' @exportMethod bind
# End generic "bind" help description in roxygen syntax
# --------------------------------------
setGeneric("bind",
           function(x,y,...){ value <- standardGeneric("bind") },
           useAsDefault=FALSE)
## --------------------------------------------
## Generic "pick" help description in roxygen syntax
## --------------------------------------------
#'  Generic function to extract one object from a complex object
#'
#' @title Extract one object from a complex object (generic function)
#'
#' @name pick
#' @aliases pick
#' @param keys an object from which to pick
#' @param selection index vector to indicate the selection
#' @return  The selected object
#' @keywords design
#' @author H. Monod, and al.
#' @seealso Method \code{\link[=pick,listofkeyrings-method]{pick-method}}
#' @exportMethod pick
# End generic "pick" help description in roxygen syntax
## --------------------------------------------------------
setGeneric("pick",
           function(keys,...){
             value <- standardGeneric("pick")
         })
## --------------------------------------------
# Generic "planor.design" help description in roxygen syntax
# --------------------------------------
#'  Generic function to build a design
#'
#' @title Build a design (generic function)
#'
#' @name planor.design
#' @aliases planor.design
#' @param key an object from which the design will be built
#' @return An object which  contains the design build from the input.
#' @keywords design
#' @author H. Monod, and al.
#' @seealso Methods \code{\link{planor.design.designkey}},
#' \code{\link{planor.design.listofkeyrings}},
#' \code{\link{planor.design.levels}}
#' @exportMethod planor.design
# End generic "design" help description in roxygen syntax
# --------------------------------------
setGeneric("planor.design",
           function(key, ...){
               value <- standardGeneric("planor.design")
           })
##---------------------------------------------------------------------------
setGeneric("getDesign",
           function(object, ...){
               value <- standardGeneric("getDesign")
           })

