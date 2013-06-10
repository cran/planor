\name{bind-methods}
\docType{methods}
\alias{bind-method}
\alias{bind}
\alias{bind,designfactors,designfactors-method}
\title{ Methods for function  \code{bind} in package \pkg{planor}: bind two objects }
\description{
Bind two objects of the same class.
}
\usage{
\S4method{bind}{designfactors,designfactors}(x, y)
}
\arguments{
  \item{x}{ an object of the first class in the signature.}
  \item{y}{ an object of the second class in the signature.}
  }
\value{An object of the same class as \code{x} and \code{y}
  containing their joint content.
}

\note{ Factors of the same name are repeated with distinct names and
  advertised with a warning.  }

\seealso{
  Class where this method applies:
  \code{\linkS4class{designfactors}}
}
\examples{
F1 <- planor.factors(factors=c("block",LETTERS[1:4]), nlevels=c(6,6,4,2,6))
F2 <- planor.factors(factors=c("block",LETTERS[11:12]), nlevels=c(6,6,4))
### Method bind on 'designfactors' objects
F3 <- bind(F1,F2)
names(F3)
}
\keyword{methods}
\keyword{ design }
