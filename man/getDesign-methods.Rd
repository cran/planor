\name{getDesign-methods}
\alias{getDesign-method}
\docType{methods}
\alias{getDesign}
\alias{getDesign,planordesign-method}
%\alias{getDesign.planordesign}

\title{  Extract a Design}
\description{
Methods to extract a design data frame from an object.
}
\usage{
\S4method{getDesign}{planordesign}( object  )
}

\arguments{
  \item{object  }{ object of the class. }
    }
\value{ A data frame
      which  contains the design.
    }
    \note{
      Applied on a \code{\linkS4class{planordesign}} object,
      it is the same as function \code{\link[=as.data.frame.planordesign]{as.data.frame}},
      without returning any attributes.
      }
      
    \seealso{
      Class where this method applies:
\code{\linkS4class{planordesign}}.
}
\examples{
### Creation of a planordesign object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
   nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
   nunits=12, base=~R+C+U, max.sol=2)
P0 <- planor.design(key=K0, select=1)
### Method getDesign on the planordesign object
show(getDesign(P0))
}
\keyword{methods}
\keyword{design}
