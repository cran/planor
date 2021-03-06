\name{listofkeyrings-class}
\Rdversion{1.1}
\docType{class}
\alias{listofkeyrings-class}
\alias{[,listofkeyrings,ANY,ANY,ANY-method}
%\alias{alias,listofkeyrings-method}
%\alias{alias.listofkeyrings}
%\alias{pick,listofkeyrings-method}
%\alias{pick.listofkeyrings}
%\alias{planor.design,listofkeyrings-method}
%\alias{planor.design.listofkeyrings}
%\alias{show,listofkeyrings-method}
%\alias{show.listofkeyrings}
%\alias{summary,listofkeyrings-method}
%\alias{summary.listofkeyrings}

\title{Class listofkeyrings}

\description{ An S4 class to store design key solutions when there is
only one prime involved or when the solutions are independent between
primes.}

\section{Objects from the Class}{ Objects are created by
 \code{\link{planor.designkey}}, when the case is not recursive.  }

\section{Slots}{
  \describe{

    \item{.Data}{a list of objects of class
      \code{\linkS4class{keyring}} associated with different primes.}

    \item{factors}{an object of class
      \code{\linkS4class{designfactors}} which contains the factors
      specifications.}

    \item{model}{a list which contains the model and
      estimate specifications.}

    \item{nunits}{the number of units of the design.}

  }
}
\section{Extends}{
Class \code{\linkS4class{list}}, from data part.
Class \code{\linkS4class{vector}}, by class list, distance 2.
}
\section{Methods}{
  \describe{
    
    \item{alias}{give the aliasing for each key-matrix. See
      \code{\link[=alias.listofkeyrings]{alias}} method.  }
    
    \item{[}{extract one design key by taking one key matrix per prime.}
    
    \item{pick}{extract one design key by taking one key matrix per
      prime.  See \code{\link{pick}} method.  }
    
    \item{planor.design}{build a design using one key matrix per prime.
      See \code{\link{planor.design}} method.  }
    
    \item{show}{display the object. See
    \code{\link[=show-method]{show}} method.  }
    
    \item{summary}{ summarize the design properties from \code{object}.
      See \code{\link[=summary.listofkeyrings]{summary}} method.

    }
  }
}

\author{Monod, H. and Bouvier, A. }

\seealso{
Creator function:  \code{\link{planor.designkey}}
}
\examples{
showClass("listofkeyrings")
### Creation of a listofkeyrings object
K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
   model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
   nunits=3^3, base=~A+B+C, max.sol=2)
show(K0)
}
\keyword{classes}
\keyword{design}
