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

\title{Class listofdesignkeys}

\description{An S4 class to represent a list of design key solutions.}

\section{Objects from the Class}{ Objects are created by
\code{\link{planor.designkey}}, when the search is recursive. }

\section{Slots}{
  \describe{
    
    \item{.Data}{a list of objects of class
      \code{\linkS4class{designkey}}.}
    
    \item{factors}{an object of class
      \code{\linkS4class{designfactors}} which contains the factors
      specifications.}
    
    \item{model}{a list which contains the model and
      estimate specifications.}
    
    \item{nunits}{the number of units in the design.}
    
  }
}

\section{Extends}{
Class \code{\linkS4class{list}}, from data part.
Class \code{\linkS4class{vector}}, by class list, distance 2.
}

\section{Methods}{
  \describe{

    \item{alias}{  NOT YET IMPLEMENTED. }
    
    \item{[}{extract one design key in the list.}
    
    \item{pick}{extract one design key in the list. See
      \code{\link[=pick]{pick}} method.}
    
    \item{planor.design}{ build a design from one design key in the
      list.  See
      \code{\link{planor.design}} method.  }
    
    \item{show}{ display the object. See
      \code{\link[=show-method]{show}} method.}
        
    \item{summary}{summarize the design properties.
      See \code{\link[=summary.listofdesignkeys]{summary}} method.
    }
  }
}

\author{Monod, H. and Bouvier, A. }
\seealso{
 Creator function:  \code{\link{planor.designkey}}
}
\examples{
showClass("listofdesignkeys")
### Creation of a listofdesignkeys object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
 nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
 nunits=12, base=~R+C+U, max.sol=2)
###  Method show
show(K0)
### Method length
length(K0)
### Extract component. The two following two commands are equivalent 
K <- K0[2]
K <- pick(K0,2)
}
\keyword{classes}
\keyword{design}