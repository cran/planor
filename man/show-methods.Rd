\name{show-methods}
\docType{methods}
\alias{show-method}
\alias{show,designkey-method}
\alias{show,keymatrix-method}
\alias{show,keyring-method}
\alias{show,listofdesignkeys-method}
\alias{show,listofkeyrings-method}

\title{  Display Objects }
\description{
 Methods to display the design key
  matrices.
}
\usage{
\S4method{show}{designkey}(object)

\S4method{show}{keymatrix}(object)

\S4method{show}{keyring}(object)

\S4method{show}{listofdesignkeys}(object)

\S4method{show}{listofkeyrings}(object)
}
\arguments{
  \item{object}{ object of the class}
}
\value{
  NULL
}
\details{
  The slot \code{pseudo.info} of the objects of class
  \code{\linkS4class{keymatrix}}  is invisible.
}

\note{
-   An R option named \code{planor.max.print} is set. It is equal
  to the number of printed
  rows and columns in the display of \pkg{planor} matrices. Default is 20.
  You can change its value by using the function
  \code{options()} (see \code{?options}).\cr
-  This method is automatically invoked when objects of the class are   displayed (see examples).
	 }
\seealso{
  Classes where this method applies:
  \code{\linkS4class{designkey}},
  \code{\linkS4class{keymatrix}},
  \code{\linkS4class{keyring}},
  \code{\linkS4class{listofdesignkeys}},
  \code{\linkS4class{listofkeyrings}}      
}
\examples{
### Creation of a listofdesignkeys object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
 nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
 nunits=12, base=~R+C+U, max.sol=2)
### Method show applied on a keymatrix object
show(K0[[1]][[1]])
### Method show applied on a designkey object
show(K0[1])
### Method show applied on the listofdesignkeys object
show(K0)
K0 #  same

### Creation of a listofkeyrings object
K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
nunits=3^3, base=~A+B+C, max.sol=2)
### Method show applied on a keyring object
show(K0[[1]]) 
print(K0[[1]]) #  same
K0[[1]] # same
### Method show applied on the listofkeyrings object
show(K0)
}

\keyword{methods}
\keyword{ design }