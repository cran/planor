\name{planordesign-class}
\Rdversion{1.1}
\docType{class}
\alias{planordesign-class}
\alias{[,planordesign,ANY,ANY,ANY-method}
%\alias{getDesign,planordesign-method}
%\alias{getDesign.planordesign}

\title{Class \code{planordesign} and methods of the class}
\description{
An S4 class to represent a final design
}
\section{Objects from the Class}{
  Objects can be created  by calls to method
\code{\link{planor.design}} applied on an object of class
  \code{\linkS4class{designkey}}
  or on an object of class
  \code{\linkS4class{listofkeyrings}}, and by calls to
    \code{\link{regular.design}} when argument \code{output} is equal to \sQuote{planordesign} }
\section{Slots}{
	 \describe{
    \item{\code{design}:}{a dataframe containing the final design}
    \item{\code{factors}:}{an object of class
      \code{\linkS4class{designfactors}} which contains the factors' specifications}
    \item{\code{model}:}{ a list containing the model and estimate
      specifications}
    \item{\code{designkey}:}{ a
	\code{list}  which contains the \code{designkey} matrices
	used to create the object}
     \item{\code{nunits}:}{the number of units of the design}
    \item{\code{recursive}:}{a \code{"logical"}
      equal to TRUE if the design has been constructed recursively}
  }
}
\section{Methods}{
  \describe{
    \item{getDesign}{ extract
      a design data frame. See  \code{\link{getDesign}}}
  }
  }
\author{
H. Monod, and A. Bouvier
}


\seealso{
  Creators: method \code{\link{planor.design}}
  applied on an object of class
  \code{\linkS4class{designkey}} or class
  \code{\linkS4class{listofkeyrings}} or class
  \code{\linkS4class{listofdesignkeys}}\cr
  See also class \code{\linkS4class{designfactors}}.
}
\examples{
showClass("planordesign")
### Creation of a 'listofdesignkeys' object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
 nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
 nunits=12, base=~R+C+U, max.sol=2)
## Creation of a 'planordesign' object from K0
P0 <- planor.design(key=K0, select=1)
show(P0)
}
\keyword{classes}
\keyword{ design }