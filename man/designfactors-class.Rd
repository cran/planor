\name{designfactors-class}
\Rdversion{1.1}
\docType{class}
\alias{designfactors-class}
\alias{[,designfactors,ANY,ANY,ANY-method}
\alias{length,designfactors-method}
\alias{names,designfactors-method}
%\alias{bind,designfactors,designfactors-method}
%\alias{bind.designfactors}

\title{Class designfactors}

\description{An S4 class to represent the design factors and  store their characteristics.}

\section{Objects from the Class}{
  Objects from this class can be created explicitly by calls to
  \code{\link{planor.factors}} or implicitly by functions such as
  \code{\link{planor.designkey}}.}

\section{Slots}{
  \describe{
    \item{fact.info}{a dataframe with one row per factor and
  with columns progressively storing information on the factors, in
  particular their numbers of levels (\code{nlev}).}
    \item{pseudo.info}{a dataframe with one row per pseudofactor
  and with columns progressively storing information on the
  pseudofactors.}
    \item{levels}{a list of numeric or character vectors, with
      each vector containing the levels of one factor.}
  }
}

\section{Methods}{
  \describe{
    \item{[}{ extract a subset of factors and update all the
  slots.}
\item{bind}{ bind two objects.
  See \code{\link[=bind]{bind}} method.}
    \item{length}{return the  number of factors. }
    \item{names}{return the names of the factors. }
	 }
}

\section{Details}{ Depending on the context and on the construction
  stage, \code{fact.info} may contain logical columns that identify the
  block factors (\code{block}), the ordered factors (\code{ordered}),
  the basic factors (\code{basic}) and so on. It may also include
  columns that store the information on the hierarchy relationships
  between factors, if any.

  In package \pkg{planor}, factors are systematically decomposed into
  pseudofactors which all have a prime number of levels and which play a
  key role in the design generation. The information on the
  pseudofactors is stored in the \code{pseudo.info} slot. In addition to
  the columns of \code{fact.info}, it contains a column (called
  \code{parent}) to give the factor that each pseudofactor
  decomposes.}

\author{ Monod, H. and Bouvier, A. }
\seealso{
 Creator function: 
  \code{\link{planor.factors}} 
}

\examples{
F1 <- planor.factors(factors=c("block",LETTERS[1:4]), nlevels=c(6,6,4,2,6))
F2 <- planor.factors(factors=c("block",LETTERS[11:12]), nlevels=c(4,6,6))
### Method bind - a warning will be issued because two factors
### in F1 and F2 have the same name
F3 <- bind(F1,F2) 
names(F3)
length(F3)
F3@levels
F3.trt <- F3[c(2:5,7,8)]
names(F3.trt)
}
\keyword{classes}
\keyword{design}
