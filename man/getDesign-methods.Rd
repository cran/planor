\name{getDesign-methods}
\alias{getDesign-method}
\alias{getDesign}
\docType{methods}
%\alias{getDesign}
\alias{getDesign,planordesign-method}
%\alias{getDesign.planordesign}

\title{ Methods for function  \code{getDesign}: extract a design}
\description{
 Extract a design data frame from an object
}
\usage{
\S4method{getDesign}{planordesign}( object  )
}

\arguments{
  \item{object  }{ Object of the class }
    }
    \value{ A data frame
      which  contains the design.
    }
    \seealso{
      Classes where this method applies:
\code{\linkS4class{planordesign}}.
}
\examples{
### Creation of a 'planordesign' object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
 nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
 nunits=12, base=~R+C+U, max.sol=2)
P0 <- planor.design(key=K0, select=1)
## Method getDesign on the 'planordesign' object
 show(getDesign(P0))
}
\keyword{methods}
\keyword{design}
