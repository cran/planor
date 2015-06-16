\name{designkey-class}
\Rdversion{1.1}
\docType{class}
\alias{designkey-class}
%\alias{alias,designkey-method}
%\alias{alias.designkey}
%\alias{planor.design,designkey-method}
%\alias{planor.design.designkey}
%\alias{show,designkey-method}
%\alias{show.designkey}
%\alias{summary,designkey-method}
%\alias{summary.designkey}

\title{Class designkey}

\description{An S4 class to represent a design-key solution. }

\section{Objects from the Class}{
Objects can be created by extraction from an object of class
\code{\linkS4class{listofkeyrings}} or class
\code{\linkS4class{listofdesignkeys}}.
}

\section{Slots}{
 \describe{
    \item{.Data}{a single design-key solution,
      i.e a list with one
      \code{\linkS4class{keymatrix}} per prime.}
    \item{factors}{an object of class
      \code{\linkS4class{designfactors}} which contains the factors specifications.}
    \item{model}{a list which  contains the model and estimate specifications.}
    \item{nunits}{the number of units of the design.}
    \item{recursive}{a logical
      equal to TRUE if the design has been constructed recursively.}
  }
}
\section{Extends}{
Class \code{\linkS4class{list}}, from data part.
Class \code{\linkS4class{vector}}, by class list, distance 2.
}
\section{Methods}{
  \describe{
    \item{alias}{ summarize the design properties.
            See \code{\link[=alias.designkey]{alias} method}.
       }
    \item{planor.design}{ 
      build the
      design from the design key matrix.
See \code{\link{planor.design} method}.
  }
  \item{show}{display the object.
    See \code{\link[=show-methods]{show} method}.
	 }
	 \item{summary}{summarize  the design properties.
	   See \code{\link[=summary.designkey]{summary} method}.
}
}}
\author{ Monod, H. and Bouvier, A. }
\examples{
### Creation of a designkey object
K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
    model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
    nunits=3^3, base=~A+B+C, max.sol=2)
print(K0[1])
}
\keyword{classes}
\keyword{design}