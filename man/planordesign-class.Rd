\name{planordesign-class}
\Rdversion{1.1}
\docType{class}
\alias{planordesign-class}
\alias{[,planordesign,ANY,ANY,ANY-method}
%\alias{getDesign,planordesign-method}
%\alias{getDesign.planordesign}

\title{Class planordesign}
\description{
An S4 class to represent a final design.
}
\section{Objects from the Class}{
  Objects can be created  by calls to method
\code{\link{planor.design}} applied on an object of class
  \code{\linkS4class{designkey}}
  or on an object of class
  \code{\linkS4class{listofkeyrings}}, and by calls to
    \code{\link{regular.design}} when argument \code{output} is equal to \sQuote{planordesign}. }
\section{Slots}{
	 \describe{
    \item{design}{a dataframe containing the final design.}
    \item{\code{factors}}{an object of class
      \code{\linkS4class{designfactors}} which contains the factors specifications.}
    \item{model}{ a list containing the model and estimate
      specifications.}
    \item{designkey}{ a
	\code{list} of the \code{\link[=keymatrix-class]{keymatrix}} objects
	used to create the object.}
     \item{nunits}{the number of units of the design.}
    \item{recursive}{a \code{logical}
      equal to TRUE if the design has been constructed recursively.}
  }
}
\section{Methods}{
  \describe{
    \item{getDesign}{ extract
      a design data frame. See  \code{\link{getDesign}} method.}
    \item{as.data.frame}{ coerce into a data frame. See
      \code{\link[=as.data.frame.planordesign]{as.data.frame}} method.}
  }
  }
\author{ Monod, H. and Bouvier, A. }



\examples{
showClass("planordesign")
### Creation of a listofdesignkeys object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
 nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
 nunits=12, base=~R+C+U, max.sol=2)
### Creation of a planordesign object from K0
P0 <- planor.design(key=K0, select=1)
show(P0)
}
\keyword{classes}
\keyword{ design }