\name{as.data.frame.planordesign}
\docType{methods}
\alias{as.data.frame.planordesign}
\alias{as.data.frame,planordesign-method}
\title{ Coerce a planordesign object to a data frame}
\description{Extracts from a \code{\linkS4class{planordesign}} object the slot \dQuote{design}, i.e the dataframe containing the final design, and stores
  the other slots in attributes}
\usage{
\S4method{as.data.frame}{planordesign}(x, \dots)
}
\value{
  A data frame with attributes \dQuote{factors}, \dQuote{model}, \dQuote{designkey}, \dQuote{nunits}, \dQuote{recursive}.}
\arguments{
  \item{x}{an object of class \code{\linkS4class{planordesign}}}
  \item{\dots}{Ignored}
}
\examples{
### Creation of a 'planordesign' object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
 nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
 nunits=12, base=~R+C+U, max.sol=2)
P0 <- planor.design(key=K0, select=1)
# Convert into a data frame
D0=as.data.frame(P0)
}

\keyword{methods}
