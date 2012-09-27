\name{keymatrix-class}
\Rdversion{1.1}
\docType{class}
\title{Class \code{'keymatrix'} and  methods of the class}
\alias{keymatrix-class}
%\alias{alias,keymatrix-method}
%\alias{alias.keymatrix}
%  \alias{show,keymatrix-method}
%  \alias{show.keymatrix}
%    \alias{summary,keymatrix-method}
%    \alias{summary.keymatrix}

\description{ An S4 class to represent an elementary key matrix in
package \pkg{planor} }

\section{Objects from the Class}{
Objects from this class are usually components of an object of class
  \code{\linkS4class{keyring}} or  \code{\linkS4class{designkey}}
}

\section{Slots}{
  \describe{
    \item{\code{.Data}:}{a matrix of integers modulo \code{p}}
    \item{\code{p}:}{a prime number}
  }
}

\section{Extends}{
Class \code{\linkS4class{matrix}}, from data part.
Class \code{\linkS4class{array}}, by class "matrix", distance 2.
Class \code{\linkS4class{structure}}, by class "matrix", distance 3.
Class \code{\linkS4class{vector}}, by class "matrix", distance 4, with explicit coerce.
}

\section{Methods}{
  \describe{
    
    \item{alias}{ gives the aliasing relationships of the key matrix.
      See \code{\link{alias-method}}.  }
    
    \item{show}{ display the object. See \link[=show-method]{show-method}}
    
    \item{summary}{ summarize the design properties.  See
      \code{\link[=summary-method]{summary-method}} }
    
  }
}

\author{
H. Monod, and A. Bouvier
}
\seealso{
 \code{\linkS4class{keyring}},  \code{\linkS4class{designkey}}
}
\examples{
showClass("keymatrix")
### Creation of a 'listofkeyrings' object
K0 <- planor.designkey(factors=c("block", LETTERS[1:4]), nlevels=rep(3,5),
  model=~block + (A+B+C+D)^2, estimate=~A+B+C+D,
  nunits=3^3, base=~A+B+C, max.sol=2)
# Method show on a 'keymatrix' of K0
show(K0[[1]][[1]])
}
\keyword{classes}
\keyword{design}
