\name{planor.model}
\title{Model and Estimate Specifications for a Design Search}
\alias{planor.model}
\usage{planor.model(model, estimate, listofmodels, resolution, factors)}

\description{Declare the factorial terms that must be
  considered as non-negligible and the factorial terms that must be
  estimable when the experiment will be analysed.}

\arguments{\item{model}{main model formula. It contains all the
    non-negligible factorial terms.}
  
  \item{estimate}{optional formula specifying the factorial terms to
    estimate. If missing, it is considered that all factorial terms in
    \code{model} have to be estimated.}
  
  \item{listofmodels}{list of \code{c(model, estimate)} pairs, where
    \code{model} and \code{estimate} are formulae; using several pairs
    allows more flexibility in the design constraints (see Kobilinsky,
    2005, or the split-plot example in the vignette); \code{estimate} is
    optional.}

  \item{resolution}{an integer larger than or equal to 3, to specify the
    design resolution. When set, the \code{model} and \code{estimate}
    arguments are ignored.  See Note.}

  \item{factors}{a \code{\linkS4class{designfactors}} object, typically
    an output from \code{\link{planor.factors}}. It must be set only when
    the \code{resolution} argument is used.}
}

\value{A list of \code{c(model, estimate)} pairs, where \code{model} and
    \code{estimate} are formulae.}

  \note{The user can specify:
    \cr 1/ either, \code{model} or \code{listofmodels} or both,
    \cr 2/ or, \code{resolution} and \code{factors}, and possibly
    \code{listofmodels}.
    
    When \code{model} and \code{resolution} are both set, \code{model} is
    ignored.

    The second case, --- when \code{resolution} and
    \code{factors} are set ---, causes the automatic generation of the
    main \code{c(model, estimate)} pair. Assuming \code{S} denotes the
    additive formula including all factors,
    \cr - if \code{resolution} is odd, the model formula is
    \code{~(S)^(resolution-1)/2},
    \cr - if \code{resolution} is even, the model formula is
    \code{~(S)^(resolution/2)} and the estimate formula is
    \code{~(S)^(resolution/2)-1}. }

\author{ Monod, H. and Bouvier, A. }
\keyword{design}
\examples{
### Basic example
planor.model(model=~block + (A+B+C)^2, estimate=~(A+B+C)^2)
### Resolution: both calls to planor.model below are equivalent
planor.model(model=~(A+B+C+D)^2, estimate=~A+B+C+D)
myfactors <- planor.factors(factors=c(LETTERS[1:4]), nlevels=rep(2,4))
planor.model(resolution=4, factors=myfactors)
### Complicated examples
planor.model(~A+B+C+D+A:B, ~A+B+C+D, listofmodels=list(c(~E+F,~E)))
planor.model(~A+B+C+D+A:B,~A+B+C+D, listofmodels=
                              list(c(~E+F,~E), ~G, ~H, c(~M+N,~N)))
}
