\name{keyring-class}
\Rdversion{1.1}
\docType{class}
\alias{keyring-class}
%\alias{show,keyring-method}
%\alias{show.keyring}
%\alias{summary,keyring-method}
%\alias{summary.keyring}
\title{Class \code{'keyring'} and methods of the class}

\description{ An S4 class to represent a list of design-key matrices
which are associated with the same prime and which represent alternative
solutions to the same design specifications.  }

\section{Objects from the Class}{
  Each component of the structure
  returned by \code{\link{planor.designkey}} is a \code{keyring} when the
  case is not recursive.  }

\section{Slots}{
  \describe{
    \item{\code{.Data}:}{a list of \code{\linkS4class{keymatrix}}
      objects.}

    \item{\code{p}:}{a prime number.}

    \item{\code{LIB}:}{a list containing a vector of row names and a
    vector of column names. The names are the same for all key matrices. }

    \item{\code{pseudo.info}:}{a dataframe containing information on the
      pseudofactors associated with the key matrices. See the
      description of the class \code{\linkS4class{designfactors}}.}

  }
}

\section{Extends}{ Class \code{"\linkS4class{list}"}, from data part.
  Class \code{"\linkS4class{vector}"}, by class "list", distance 2.  }

\section{Methods}{ \describe{

    \item{show}{ display the object. See \link[=show-method]{show-method}.}

    \item{summary}{ summarize the design properties. See
      \code{\link[=summary-method]{summary-method}}}

  }
}

\author{
 H. Monod, and A. Bouvier
}

\note{Each key matrix in a \code{keyring} object is a possible solution
  to the same factors, model and estimate specifications, with respect
  to the same prime number. An object of class
  \code{\linkS4class{listofkeyrings}} is a list of \code{keyring}
  objects associated with the different primes involved in a given
  factorial design problem.  }

\seealso{
\code{\link{planor.designkey}},
method \code{pick} in class \code{\linkS4class{listofkeyrings}},
method \code{summary} in class
\code{\linkS4class{keymatrix}} and the class \code{\linkS4class{keyring}}
}
\examples{
showClass("keyring")
### Creation of a 'listofkeyrings' object
K0 <- planor.designkey(factors=c("block", LETTERS[1:4]), nlevels=rep(3,5),
  model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
  nunits=3^3, base=~A+B+C, max.sol=2)
## Method show applied on a 'keyring' component of K0
show(K0[[1]]) 
}
\keyword{classes}
\keyword{design}