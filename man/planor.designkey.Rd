\name{planor.designkey}
\title{Search for a Design Key or a Collection of Design Keys}
\alias{planor.designkey}

\description{Search for a design key or a collection of design keys that
satisfy the design properties specified by the arguments. This function
implements the core algorithms of the  \pkg{planor} package.}

\usage{planor.designkey(factors, nlevels, block, ordered, hierarchy, model,
    estimate, listofmodels, resolution, nunits, base, max.sol=1,
    randomsearch=FALSE, verbose=TRUE)}

\arguments{
  \item{factors}{an object of class \code{\linkS4class{designfactors}},
    typically an output from
  \code{\link{planor.factors}}. Alternatively, you can use
    the arguments \code{factors}, \code{nlevels}, \code{ordered},
    \code{hierarchy} as described in the syntax of \code{\link{planor.factors}}.}
  
  \item{nlevels}{see \code{\link{planor.factors}}. Ignored if
    \code{factors} is of class \code{\linkS4class{designfactors}}.}

  \item{block}{see \code{\link{planor.factors}}. Ignored if
    \code{factors} is of class \code{\linkS4class{designfactors}}.}
  
  \item{ordered}{see \code{\link{planor.factors}}. Ignored if
    \code{factors} is of class \code{\linkS4class{designfactors}}.}
  
  \item{hierarchy}{see \code{\link{planor.factors}}. Ignored if
    \code{factors} is of class \code{\linkS4class{designfactors}}.}
  
  \item{model}{a list of model-estimate pairs of formulae, typically an
    output from \code{\link{planor.model}}. Alternatively, you can use
    the arguments
    \code{model}, \code{estimate}, \code{listofmodels} and
    \code{resolution}, as described in the syntax of
    \code{\link{planor.model}}.}
  
  \item{estimate}{see \code{\link{planor.model}}. Ignored if
    \code{model} is a list.}
  
  \item{listofmodels}{see \code{\link{planor.model}}. Ignored if
    \code{model} is a list.}

  \item{resolution}{see \code{\link{planor.model}}. Ignored if
    \code{model} is a list. When  set and there is no hierarchy, a faster algorithm is used
    which exploits the  symmetries. }
  
  \item{nunits}{a scalar giving the total number of units in the design}
  
  \item{base}{an optional additive formula to specify the basic
    factors. See Note.}
  
  \item{max.sol}{maximum number of solutions before exit.}
  
  \item{randomsearch}{a \code{logical}. If TRUE, the searches for a key
  matrix are performed in a random order.}
  
  \item{verbose}{a \code{logical} to set to \code{TRUE} for verbose
  display.}}

\value{An object of class \code{\linkS4class{listofkeyrings}} in most
  cases. Otherwise, i.e in recursive cases, an object of class
  \code{\linkS4class{listofdesignkeys}}.}

\note{The \code{nunits} argument is compulsory except if the \code{base}
  argument is used. When both arguments are missing, the program stops and
  it gives the size that would be required by a full factorial
  design. When only \code{nunits} is missing, the number of units is
  given by the product of the numbers of levels of the base
  factors.

  The \code{base} formula must be an additive formula involving a
  subset of factors, called the basic factors. Using the \code{base}
  argument ensures that the design solutions will include the full
  factorial design for the basic factors. This option can speed up
  the search because it restricts the domain to be explored by the search
  algorithm. }

\details{The methods implemented in \pkg{planor} rely on a decomposition
  of the design search according to prime numbers. The prime numbers
  involved are those that decompose the numbers of levels of the
  factors. For example, if all factors have 2, 4, or 8 levels, then the
  number of units must be a power of 2 and the only prime number
  involved is 2. This is called the \emph{symmetric} case. But if at
  least one factor has 6 levels, or if factor \eqn{A} has 2 levels and
  factor \eqn{B} has 3 levels, then the number of units must be the
  product of a power of 2 by a power of 3. In this case the search is
  automatically decomposed into one for prime 2 and one for prime
  3. This is called the \emph{asymmetric} case.

  In the symmetric case with prime \eqn{p}, a regular factorial design
  requires a single key matrix of integers modulo \eqn{p}. In the
  asymmetric case, it requires one key matrix per prime. In
  \pkg{planor}, key matrices are stored in objects of class
  \code{\linkS4class{keymatrix}}. The lists made of one key matrix per
  prime are called design keys. They are stored in objects of class
  \code{\linkS4class{designkey}}.

  The function \code{\link{planor.designkey}} essentially searches for
  design keys that satisfy the user specifications. For technical
  reasons, however, its output can take two different forms: either an
  object of class \code{\linkS4class{listofkeyrings}} or an object of
  class \code{\linkS4class{listofdesignkeys}}. The function
  \code{\link{planor.designkey}} detects automatically which case
  applies. In the first case (\emph{independent case}), the key matrix
  solutions can be searched independently between primes and they are
  stored in objects of class \code{\linkS4class{listofkeyrings}}. The
  second case (\emph{recursive case}) occurs exceptionnally. In that
  case the search cannot be independent between primes and so the
  different solutions are directly stored in a list of class
  \code{\linkS4class{listofdesignkeys}}.  }

\author{ Monod, H. and Bouvier, A. }

\seealso{\code{\link{planor.factors}}, \code{\link{planor.model}}, and
  the classes \code{\linkS4class{designfactors}},
  \code{\linkS4class{listofkeyrings}},
  \code{\linkS4class{listofdesignkeys}}}
\keyword{design}
\examples{
K0 <- planor.designkey(factors=c("block", LETTERS[1:4]),
  nlevels=rep(3,5), model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
  nunits=3^3, base=~A+B+C, max.sol=2)
### With automatic model generation
Km <- planor.designkey(factors=c("block", LETTERS[1:4]),
  nlevels=rep(2,5), resolution=3, nunits=2^4)
}
