\name{bind-methods}
\docType{methods}
\alias{bind,designfactors,designfactors-method}
\alias{bind}
\title{  Bind Two Objects }
\description{
Methods to bind two objects of the same class.
}
\usage{
\S4method{bind}{designfactors,designfactors}(x, y)
}
\arguments{
  \item{x}{ an object of the first class in the signature.}
  \item{y}{ an object of the second class in the signature.}
  }
\value{An object of the same class as \code{x} and \code{y},
  containing their joint content.
}

\note{ Factors with same name are  distinguished and
  advertised with a warning.  }

\seealso{
  Class where this method applies:
  \code{\linkS4class{designfactors}}
}
\examples{
F1 <- planor.factors(factors=c("block",LETTERS[1:4]), nlevels=c(6,6,4,2,6))
F2 <- planor.factors(factors=c("block",LETTERS[11:12]), nlevels=c(6,6,4))
### Method bind on designfactors objects
F3 <- bind(F1,F2)
names(F3)
}
\keyword{methods}
\keyword{ design }
