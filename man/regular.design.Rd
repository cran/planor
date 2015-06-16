\name{regular.design}
\title{Construct and Randomize a Regular Factorial Design}
\usage{
regular.design(factors = NULL, nlevels = NULL, block = NULL,
   ordered = NULL, hierarchy = NULL, model = NULL, estimate = NULL,
   listofmodels = NULL, resolution = NULL, nunits = NULL,
   base = NULL, randomize = NULL, randomsearch = FALSE,
   output = "planordesign", verbose = FALSE, ...)
}
\description{Construct and randomize a regular factorial design.}
\alias{regular.design}
\arguments{
  \item{factors}{an object of class \code{\linkS4class{designfactors}},
    typically an output from \code{\link{planor.factors}}). Otherwise
    the arguments \code{factors}, \code{nlevels}, \code{ordered},
    \code{hierarchy} follow the syntax of \code{\link{planor.factors}}.}
  
  \item{nlevels}{see \code{\link{planor.factors}}.  Ignored if
    \code{factors} is of class \code{\linkS4class{designfactors}}.}

  \item{block}{see \code{\link{planor.factors}}.  Ignored if
    \code{factors} is of class \code{\linkS4class{designfactors}}.}
  
  \item{ordered}{see \code{\link{planor.factors}}.  Ignored if
    \code{factors} is of class \code{\linkS4class{designfactors}}.}
  
  \item{hierarchy}{see \code{\link{planor.factors}}.  Ignored if
    \code{factors} is of class \code{\linkS4class{designfactors}}.}
  
  \item{model}{a list of model-estimate pairs of formulae, typically an
    output from \code{\link{planor.model}}. Otherwise the arguments
    \code{model}, \code{estimate}, \code{listofmodels} and
    \code{resolution} follow the syntax of
    \code{\link{planor.model}}.}
  
  \item{estimate}{see \code{\link{planor.model}}. Ignored if
    \code{model} is a list.}
  
  \item{listofmodels}{see \code{\link{planor.model}}. Ignored if
    \code{model} is a list.}

  \item{resolution}{see \code{\link{planor.model}}. Ignored if
    \code{model} is a list.}
  
  \item{nunits}{see \code{\link{planor.designkey}}.}
  
  \item{base}{see \code{\link{planor.designkey}}.}
  
  \item{randomize}{an optional formula to randomize the design.}

  \item{randomsearch}{see \code{\link{planor.designkey}}.}

  \item{output}{a string to specify the class of the output value: either a
    \code{data.frame} or a \code{\link[=planordesign-class]{planordesign}} object.}

  \item{verbose}{a \code{logical} to set to \code{TRUE} for verbose
  display.}

  \item{...}{additional arguments, in particular those related to
    randomization}
}

\value{An object of class \code{data.frame} or
  \code{\linkS4class{planordesign}}, depending on the \code{output}
  argument.}

\author{Monod, H. and Bouvier, A. }
\seealso{\code{\link{planor.factors}}, \code{\link{planor.model}},
  and the classes  \code{\linkS4class{designfactors}},
  \code{\linkS4class{listofkeyrings}},
  \code{\linkS4class{listofdesignkeys}}}
\keyword{design}
\examples{mydesign <- regular.design(factors=c("block", LETTERS[1:4]),
  nlevels=rep(3,5), model=~block + (A+B+C+D)^2, estimate=~A+B+C+D,
  nunits=3^3, randomize=~block/UNITS)
print(mydesign)
}
