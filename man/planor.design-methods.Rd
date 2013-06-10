\name{planor.design-methods}
\docType{methods}
\alias{planor.design}
\alias{planor.design,designkey-method}
\alias{planor.design,listofdesignkeys-method}
\alias{planor.design,listofkeyrings-method}
\alias{planor.design,numeric-method}

\title{ Methods for function \code{planor.design} in package  \pkg{planor}: build a design from
a design-key solution}

\description{Construction of a factorial design from an object
  containing key matrices}

\usage{
\S4method{planor.design}{designkey}(key, randomize=NULL, \dots)

\S4method{planor.design}{listofdesignkeys}(key, randomize=NULL, selection=1, \dots)

\S4method{planor.design}{listofkeyrings}(key, randomize=NULL, selection,\dots)

\S4method{planor.design}{numeric}(key, start=1)

}

\arguments{

  \item{key}{an object of the first class in the signature, or a vector
  of integers.}

  \item{randomize}{an optional formula to specify the block structure
    for design randomization.}
  
  \item{selection}{
    when \code{key} is a \code{\linkS4class{listofdesignkeys}} object,
    an integer scalar.\cr
    when \code{key} is a \code{\linkS4class{listofkeyrings}} object,
    should be an index vector to select the key matrix for each prime.
  }
  
  \item{\ldots}{additional arguments, in particular those related to
    randomization (see \code{\link{planor.randomize}}).}
  
  \item{start}{ an integer from where to start the series of symbols.}
  
}


\value{An object of class \code{\linkS4class{planordesign}}, which
  contains the design built from the input. This function is restricted
  to giving a single design.  When \code{key} is numeric, see Details}


\details{ \itemize{

    \item When \code{key} is numeric, it should be a vector of integers
    of length \code{s}.  Then, the function generates a full factorial \eqn{n_1
    x n_2 x ... x n_s} design with columns considered as factors. It returns
    an integer matrix with prod(\eqn{n}) rows and \eqn{s} columns giving all
    combinations along the rows, in lexicographic order.
    
    \item When \code{key} is a \code{\linkS4class{listofdesignkeys}}
    object, build one design from a selected solution.  } }

\seealso{
  Classes where this method applies:
  \code{\linkS4class{designkey}},
  \code{\linkS4class{listofdesignkeys}},
  \code{\linkS4class{listofkeyrings}}.
}

\examples{
### Creation of a 'listofdesignkeys' object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
 nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
 nunits=12, base=~R+C+U, max.sol=2)
## Method planor.design applied on the 'listofdesignkeys' object
P0 <- planor.design(key=K0, select=1)
## Method planor.design applied on a designkey' object
P0 <- planor.design(K0[1])


### Creation of a 'listofkeyrings' object
K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
nunits=3^3, base=~A+B+C, max.sol=2, verbose=TRUE)
## Method planor.design applied on a designkey' object
P0 <- planor.design(K0[1])
P0.R <- planor.design(K0[1], randomize=~A+B+C+D) ## randomize the final design}

\keyword{methods}
\keyword{design}
