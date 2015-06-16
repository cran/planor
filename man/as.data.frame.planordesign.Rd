\name{as.data.frame.planordesign}
\docType{methods}
\alias{as.data.frame.planordesign}
\alias{as.data.frame,planordesign-method}
\title{  Data Frame Coercion}
\description{Method to extract  the dataframe
  containing
  the final design from a
  \code{\linkS4class{planordesign}} object,
  i.e the slot \code{design}.
  The other slots are stored in attributes.}
\usage{
\S4method{as.data.frame}{planordesign}(x, \dots)
}
\value{
  A data frame with attributes \code{factors}, \code{model}, \code{designkey}, \code{nunits}, \code{recursive}.}
\arguments{
  \item{x}{an object of class \code{\linkS4class{planordesign}}.}
  \item{\dots}{ignored.}
}
\examples{
### Creation of a planordesign object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
 nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
 nunits=12, base=~R+C+U, max.sol=2)
P0 <- planor.design(key=K0, select=1)
### Convert into a data frame
D0 <- as.data.frame(P0)
}

\keyword{methods}
