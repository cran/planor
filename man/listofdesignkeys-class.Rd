\name{listofdesignkeys-class}
\Rdversion{1.1}
\docType{class}
\alias{listofdesignkeys-class}
\alias{[,listofdesignkeys,ANY,ANY,ANY-method}
%\alias{alias,listofdesignkeys-method}
%\alias{alias.listofdesignkeys}
%  \alias{planor.design,listofdesignkeys-method}
%\alias{planor.design.listofdesignkeys}
 %   \alias{show,listofdesignkeys-method}
%\alias{show.listofdesignkeys}
  %    \alias{pick,listofdesignkeys-method}
%\alias{pick.listofdesignkeys}
  %    \alias{summary,listofdesignkeys-method}
%\alias{summary.listofdesignkeys}

\title{Class \code{listofdesignkeys} and methods of the class}

\description{An S4 class to represent a list of design key solutions}

\section{Objects from the Class}{ Objects are created by
\code{\link{planor.designkey}}, when the search is recursive. }

\section{Slots}{
  \describe{
    
    \item{\code{.Data}:}{a list of objects of class
      \code{\linkS4class{designkey}}.}
    
    \item{\code{factors}:}{an object of class
      \code{\linkS4class{designfactors}} which contains the factors'
      specifications.}
    
    \item{\code{model}:}{a \code{"list"} which contains the model and
      estimate specifications.}
    
    \item{\code{nunits}:}{the number of units in the design.}
    
  }
}

\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}

\section{Methods}{
  \describe{

    \item{alias}{ FUNCTION NOT YET IMPLEMENTED }
    
    \item{[}{extract one design key in the list.}
    
    \item{pick}{extract one design key in the list. See
      \code{\link[=pick]{pick}}.}
    
    \item{planor.design}{ build a design from one design key in the
      list.  See
      \code{\link{planor.design}}.  }
    
    \item{show}{ display the object. See
      \code{\link[=show-method]{show}}.}
        
    \item{summary}{summarize the design properties.
      See \code{\link[=summary.listofdesignkeys]{summary}}
    }
  }
}

\author{
H. Monod, and A. Bouvier
}
\seealso{
 Creator function:  \code{\link{planor.designkey}}
}
\examples{
showClass("listofdesignkeys")
### Creation of a "listofdesignkeys" object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
 nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
 nunits=12, base=~R+C+U, max.sol=2)
# Show the object
show(K0)
## Method length
length(K0)
## Extraction: the following two commands are equivalent 
K <- K0[2]
K <- pick(K0,2)
}
\keyword{classes}
\keyword{design}