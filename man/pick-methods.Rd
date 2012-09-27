\name{pick-methods}
\docType{methods}
\alias{pick-method}
\alias{pick}
%\alias{pick.listofkeyrings}
\alias{pick,listofkeyrings-method}
%\alias{pick.listofdesignkeys}
\alias{pick,listofdesignkeys-method}
\title{ Methods for function  \code{pick}: extract a single  \code{designkey} from a list }
\description{
Extract a single \code{\linkS4class{designkey}}
object (with one key matrix per prime) from a complex object
}
\usage{
\S4method{pick}{listofdesignkeys}(keys, selection)

\S4method{pick}{listofkeyrings}(keys, selection)
}
\arguments{
  \item{keys}{ an object of the class }
  \item{selection}{
     when \code{keys} is a \code{\linkS4class{listofdesignkeys}} object,
     an integer scalar equal to the position of the
    required solution.\cr
    when \code{keys} is a \code{\linkS4class{listofkeyrings}} object,
    the index vector to select the key matrix for each prime.
   }
 }
  
\value{An object of class \code{\linkS4class{designkey}}, which contains
  the selected design}

\seealso{
  Classes where this method applies:
  \code{\linkS4class{listofdesignkeys}},
  \code{\linkS4class{listofkeyrings}}.
}


\note{\code{K <- pick(K0,1)} can be simply written \code{K <- K0[1]}
}
\examples{
# Creation of an object of class "listofdesignkeys"
K2 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
nlevels=c(3,2,2,3,2,2),  model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2 , nunits=12,
base=~R+C+U, max.sol=2)
# Method 'pick' applied on the "listofdesignkeys" object
K2.1 <- pick(K2,1)
K2.1 <- K2[1] ## Another way of extracting ([ is synonym of pick)

# Creation of an object of class "listofkeyrings"
K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"),
nlevels=rep(3,5), model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
nunits=3^3, base=~A+B+C, max.sol=2)
# Method 'pick' applied on the "listofkeyrings" object
K0.1 <- pick(K0,1)
K0.1 <- K0[1] ## the same
}

\keyword{methods}
\keyword{ design}
