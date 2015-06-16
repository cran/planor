\name{alias-methods}
\docType{methods}
\alias{alias,designkey-method}
\alias{alias.designkey}
\alias{alias,keymatrix-method}
\alias{alias.keymatrix}
\alias{alias,listofdesignkeys-method}
\alias{alias.listofdesignkeys}
\alias{alias,listofkeyrings-method}
\alias{alias.listofkeyrings}
\alias{alias,planordesign-method}
\alias{alias.planordesign}
\title{   Summarize the Design Properties  }
\description{
  Methods to summarize the design properties of an object containing key matrices.
  Display the design keys matrix(ces) and the factorial effects confounded with the mean.
}
\usage{
\S4method{alias}{designkey}(object, model, \dots)

\S4method{alias}{keymatrix}(object, model, fact, block, \dots)

\S4method{alias}{listofdesignkeys}(object, model, \dots)

\S4method{alias}{listofkeyrings}(object, model, \dots)

\S4method{alias}{planordesign}(object, model, fact, block, \dots)

}
\arguments{
  \item{object}{an object of the class.}
  \item{model}{an optional model formula (by default the first model in
    object) or, when \code{object} is a \code{\linkS4class{keymatrix}},
      a matrix representing factorial model terms.}
  \item{fact}{	a character or numeric vector of parent factor names for
    the columns of \code{object}.}
  \item{block}{	a logical vector to identify the columns of \code{object}
    associated with a block factor.}
  \item{ \ldots }{     ignored.}
}
  \value{
    \itemize{
    \item When \code{object} is a \code{\linkS4class{keymatrix}},
    a vector with
(i) the number of unaliased treatment effecs; (ii) the number of
mutually aliased treatment effects; (iii) the number of treatment effects
aliased with block effects.
\item When \code{object} is a \code{\linkS4class{designkey}},
 an invisible NULL.
 \item When \code{object} is a \code{\linkS4class{listofkeyrings}},
 the factors, the model and the number of solutions for each prime
in a  list indexed by the primes p of the object. Each element is a 3-column
matrix with one row per solution for prime p. The columns
give (i) the number of unaliased treatment effecs; (ii) the number of
mutually aliased treatment effects;  (iii) the number of treatment effects
aliased with block effects.
\item The method is NOT YET IMPLEMENTED on objects of class
\code{\linkS4class{listofdesignkeys}}.
\item When \code{object} is a \code{\linkS4class{planordesign}},
this function is the alias method applied on each of the \code{\linkS4class{keymatrix}}
objects contained in its \code{designkey} slot.

 }
}
\details{
  When \code{object} is a \code{\linkS4class{keymatrix}},
    \dQuote{alias} displays the key matrix and the factorial effects confounded with
      the mean. It prints the unaliased treatment effects,
then the groups of aliased treatment effects, then the treatments effects confounded with block effects
and finally the unaliased block effects, when considering all the factorial terms that are represented in the \code{model} argument,
which is set if missing to the identity matrix (main effects only).
}
\examples{
### Creation of an object of class listofkeyrings
K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
nunits=3^3, base=~A+B+C, max.sol=2)
### alias on an object of class keymatrix
alias(K0[[1]][[1]])
### alias on an object of class designkey
alias(K0[1])
### alias on an object of class listofkeyrings
alias(K0)
}


\keyword{methods}
\keyword{design}
