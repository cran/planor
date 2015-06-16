\name{planor.harmonize}
\title{Harmonize the Factors}
\alias{planor.harmonize}
\usage{planor.harmonize(factors, nlevels, ordered, hierarchy, model, estimate,
    listofmodels, base)}

\description{Harmonize the factors originating from a list of factors, a
  list of models, and a list of basic factors (this function is
  essentially for internal use).}

\arguments{
  \item{factors}{an object of class \code{\linkS4class{designfactors}},
    typically an output from \code{\link{planor.factors}}). Otherwise
    the arguments \code{factors}, \code{nlevels}, \code{ordered},
    \code{hierarchy} follow the syntax of \code{\link{planor.factors}}.}
  
  \item{nlevels}{see \code{\link{planor.factors}}.  Ignored if
    \code{factors} is of class \code{\linkS4class{designfactors}}.}

  \item{ordered}{see \code{\link{planor.factors}}. Ignored if
    \code{factors} is of class \code{\linkS4class{designfactors}}.}

  \item{hierarchy}{see \code{\link{planor.factors}}. Ignored if
    \code{factors} is an object of class
    \code{\linkS4class{designfactors}}.}

  \item{model}{a list of model-estimate pairs of formulae, typically an
    output from \code{\link{planor.model}}. Otherwise the arguments
    \code{model}, \code{estimate} and \code{listofmodels} follow
    the syntax of \code{\link{planor.model}}.}

  \item{estimate}{see \code{\link{planor.model}}. Ignored if
    \code{model} is a list.}

  \item{listofmodels}{see \code{\link{planor.model}}. Ignored if
    \code{model} is a list.}

  \item{base}{an optional formula to specify the basic factors. These
  factors must belong to the \code{factors} argument}}

\value{An object of class \code{\linkS4class{designfactors}} very
  similar to \code{factors}, but with two additional logical columns in
  slots \code{fact.info} and \code{pseudo.info}:\cr
- \code{model} (TRUE for the factors present in at least one model formula),\cr
- \code{basic} (TRUE for the basic factors).}

\note{This function is called at the start of the design search. It is
essentially a check that the factors in all three arguments are
coherent, even though it performs some additional tasks.  The function
stops if it detects a model or basic factor that is absent from
\code{factors}. This is because the number of levels of such a factor is
unknown and so the design search cannot proceed.  Besides, the function
eliminates the factors that do appear neither in \code{model} nor in
\code{base} and it reorders the factors by putting first the basic
ones.}

\keyword{design}
\examples{
F2 <- planor.factors(factors=c("block",LETTERS[1:4]), nlevels=c(6,6,6,4,2))
M2 <- planor.model( model=~block+(A+B+C)^2, estimate=~A+B+C )
F2.h <- planor.harmonize(factors=F2, model=M2, base=~A+B)
names(F2)
names(F2.h)
}
\author{ Monod, H. and Bouvier, A. }
\keyword{design}
