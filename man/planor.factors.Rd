\name{planor.factors}
\title{Create an object of class 'designfactors'}
\alias{planor.factors}
\usage{planor.factors(factors = NULL, nlevels = NULL,
  block = NULL, ordered = NULL, hierarchy = NULL,
  dummy = FALSE)}
\description{A function to create an
  object of class \code{\linkS4class{designfactors}},
  either by giving the factor names and their numbers of levels, or by giving
  a named list of factor levels. Both ways can be used in the same
  call. Additional information can be provided that will be used during the
  design search or in the summary functions.}
\arguments{\item{factors}{a character vector of factor names, or
    possibly a scalar, a dataframe or a list (see DETAILS)}
  
  \item{nlevels}{a vector of level numbers for each factor name (see
    DETAILS)}
  
  \item{block}{an additive model formula to indicate the block factors}
  
  \item{ordered}{an additive model formula to indicate the quantitative
    factors (not used at all in the present version)}
  
  \item{hierarchy}{a formula or a list of formulae to indicate hierarchy
    relationships between factors (see the \pkg{planor} vignette for
    details)}
  
  \item{dummy}{a logical to identify dummy factors created and deleted by PLANOR
    functions for technical reasons}}

\value{An object of class \code{\linkS4class{designfactors}}}

\note{The basic usage is to specify the names of the factors by a
  character vector of length \eqn{n} in argument \code{factors} and
  their numbers of levels by a numeric vector of length \eqn{n} in
  argument \code{nlevels}. Alternatively, the \code{factors} argument
  can be an integer \eqn{n}, in which case the first \eqn{n} capital
  letters of the alphabet are used as factor names.  If \code{nlevels}
  is a scalar \eqn{s}, it is considered that all factors have \eqn{s}
  levels.  There are two more possibilities which allow for alphanumeric
  factor levels.  If \code{factors} is a dataframe, the factors in this
  dataframe are extracted together with their levels. Finally
  \code{factors} can be a named list of \eqn{n} vectors, with each
  vector containing the levels of one factor. Note that \code{nlevels}
  is ignored in these latter two cases. See the examples. The argument
  \code{block} allows to specify the block or nuisance factors. This
  information is used by the \code{\link{alias}} and
  \code{summary} functions but it has no effect on the design
  generation and randomization which depend on other arguments.}

\author{H. Monod, and A. Bouvier}
\seealso{Class \code{\linkS4class{designfactors}}}
\keyword{design}
\examples{planor.factors(c("A","B","C","P"),c(2,3,6,3))
planor.factors(LETTERS[1:12],2)
planor.factors(12,2)
planor.factors( c("A","B","Block"), 3, block=~Block )
zz <- planor.factors( c("A","B","Block"), c(2,3,5))
zz@levels$A <- c("plus","moins")
planor.factors(factors=list(A=c("plus","moins"), B=1:3, Block=1:5))
AB <- data.frame( A=c(rep(c("a","b"),3)), B=rep(c("z","zz","zzz"),rep(2,3)), C=1:6  )
planor.factors(factors=AB)}
