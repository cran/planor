\name{keymatrix-class}
\Rdversion{1.1}
\docType{class}
\alias{keymatrix-class}
%\alias{alias,keymatrix-method}
%\alias{alias.keymatrix}
%  \alias{show,keymatrix-method}
%  \alias{show.keymatrix}
%    \alias{summary,keymatrix-method}
%    \alias{summary.keymatrix}

\title{Class keymatrix}
\description{ An S4 class to represent an elementary key matrix. }

\section{Objects from the Class}{
Objects from this class are usually components of an object of class
  \code{\linkS4class{keyring}} or  \code{\linkS4class{designkey}.}
}

\section{Slots}{
  \describe{
    \item{.Data}{a matrix of integers modulo \code{p}.}
    \item{p}{a prime number.}
  }
}

\section{Extends}{
Class \code{\linkS4class{matrix}}, from data part.
Class \code{\linkS4class{array}}, by class matrix, distance 2.
Class \code{\linkS4class{structure}}, by class matrix, distance 3.
Class \code{\linkS4class{vector}}, by class matrix, distance 4, with explicit coerce.
}

\section{Methods}{
  \describe{
    
    \item{alias}{ give the aliasing relationships of the key matrix.
      See \code{\link[=alias.keymatrix]{alias}} methods.  }
    
    \item{show}{ display the object.
      See \code{\link[=show-method]{show}} methods.  }
    
    \item{summary}{ summarize the design properties.
      See \code{\link[=summary.keymatrix]{summary}} methods. }
    
  }
}

\author{Monod, H. and Bouvier, A. }

\seealso{
 \code{\linkS4class{keyring}},  \code{\linkS4class{designkey}}
}
\examples{
showClass("keymatrix")
### Creation of a listofkeyrings object
K0 <- planor.designkey(factors=c("block", LETTERS[1:4]), nlevels=rep(3,5),
    model=~block + (A+B+C+D)^2, estimate=~A+B+C+D,
    nunits=3^3, base=~A+B+C, max.sol=2)
### Show a keymatrix of K0
show(K0[[1]][[1]])
}
\keyword{classes}
\keyword{design}
