\name{summary-methods}
\docType{methods}
%\alias{summary-method}
\alias{summary,designkey-method}
\alias{summary.designkey}
\alias{summary,keymatrix-method}
\alias{summary.keymatrix}
\alias{summary,keyring-method}
\alias{summary.keyring}
\alias{summary,listofdesignkeys-method}
\alias{summary.listofdesignkeys}
\alias{summary,listofkeyrings-method}
\alias{summary.listofkeyrings}
\alias{summary,planordesign-method}
\alias{summary.planordesign}

\title{  Summarize the Design Properties }
\description{
 Methods to summarize the design properties of an object, by
printing the summary of each key matrix.
}
\usage{
\S4method{summary}{designkey}(object, show="dtbw", save="k", \dots)

\S4method{summary}{keymatrix}(object, fact, block,  show="dtbw", save="k", \dots)

\S4method{summary}{keyring}(object, show="tbw",save ="kw", \dots)

\S4method{summary}{listofdesignkeys}(object, show="tbw", save="kw", \dots)

\S4method{summary}{listofkeyrings}(object,  show="tbw",  save="kw", \dots)

\S4method{summary}{planordesign}(object, fact, block,  show="dtbw", save="k", \dots)
}
\note{
  An R option named \code{planor.max.print} is set. It is equal
  to the number of printed
  rows and columns in the display of \pkg{planor} matrices. Default is 20.
  You can change its value by using the function
  \code{options()} (see \code{?options}).
  }
\value{
  A list. See Details\cr
  Information returned for each key matrix
  depends on the argument \code{save}.
  \itemize{
    \item When \code{save} includes the character \sQuote{k},
    the returned list has a component named \sQuote{k}.It is a matrix,
    the columns of which are kernel generators of the key matrices.
   \item When \code{save} includes the character \sQuote{w},
the returned list has a component named \sQuote{w}, which contains
the  weight profiles of the
	effects confounded with the mean.
}
}


\details{
  The amount of display depends on the   value
  of the argument \code{show}, and the type of returned information
    depends on the   value
   of the argument \code{save}.
   
  \itemize{
 \item When \code{object} is of class
    \code{\linkS4class{keymatrix}},
    \dQuote{summary} prints the key matrix, the factorial
	effects confounded with the mean, and the weight profiles of the
	effects confounded with the mean, according to the value
	of the argument \code{show}.\cr
	The \code{\linkS4class{keymatrix}} argument being denoted by \code{key},\cr
- The rows of \code{key} are associated with units factors (or pseudofactors)
while its columns are associated with treatment or block factors (or pseudofactors).\cr
- The vectors in the arguments \code{fact} and \code{block} give information on the treatment and block factors,
so that their length is expected to be equal to the number of columns of \code{key}.\cr
- If missing, \code{fact} attributes a distinct parent factor to each
column of \code{key} and \code{block} is set to TRUE for all columns.

\dQuote{summary} returns a list with the components required by the  argument \code{save}.
\item When \code{object} is of class \code{\linkS4class{designkey}},
\dQuote{summary} prints the summary of each of the key matrices.
It returns a list with as many components as key matrices, each one with
the components required by the  argument \code{save}.

    \item When \code{object} is of class
    \code{\linkS4class{listofdesignkeys}},
    \dQuote{summary} prints the summary of each key matrix in each
    design key. It returns a list with as many components as design keys,
    each one is a list of the  key matrices summaries. 
    \item When \code{object} is of class
    \code{\linkS4class{listofkeyrings}},
    \dQuote{summary} prints  the summary of each key matrix in each
    keyring.  It returns a list with as many components as keyrings,
    each one is a list of the  key matrices summaries. 
    \item When \code{object} is of class
    \code{\linkS4class{keyring}},
    \dQuote{summary} prints   the summary of each of its key matrices.
    It returns a list with as many components as  key matrices.
\item When \code{object} is a \code{\linkS4class{planordesign}},
this function is the summary method applied on each of the  keymatrix
objects contained in its \code{designkey} slot.

	
  }
}

\arguments{
  \item{object}{an object of the class.}
  \item{fact}{a character or numeric vector of parent factor names for
    the columns of the object.}
\item{block}{a logical vector to identify the columns of the object
  associated with a block factor}
\item{show}{ an optional string to  identify the
  type of information to display. The recognized letters are: \sQuote{d} for the design keys matrices, \sQuote{t}
  for the treatment effects
  confounded with the mean,
  \sQuote{b} for the block-and-treatment effects confounded with the mean,
  \sQuote{w} for the weight profiles.}
\item{save}{ an optional string to  identify the
  type of information to return. The recognized letters are:
  \sQuote{k} for the kernel matrices, \sQuote{w} for the  weight profiles of the
treatment effects confounded with the mean.}
\item{\ldots}{ignored.}
}
\examples{
### Creation of a listofdesignkeys object
K0 <- planor.designkey(factors=c("R","C","U","A","B1","B2"),
  nlevels=c(3,2,2,3,2,2), model=~R*C + (A+B1+B2)^2, estimate=~A:B1+A:B2,
  nunits=12, base=~R+C+U, max.sol=2)
### Method summary applied on a keymatrix object
r <- summary(K0[[1]][[1]])
### Method summary applied on a designkey  object
summary(K0[1], save=NULL)
### Method summary applied on the listofdesignkeys object
r <-summary(K0, show="dt")

### Creation of a listofkeyrings object
K0 <- planor.designkey(factors=c(LETTERS[1:4], "block"), nlevels=rep(3,5),
   model=~block+(A+B+C+D)^2, estimate=~A+B+C+D,
   nunits=3^3, base=~A+B+C, max.sol=2)
### Method summary applied on the keymatrix object
r <-summary(K0[[1]][[1]])
### Method summary applied on the keyring object
r <-summary(K0[[1]])
### Method summary applied on the listofkeyrings object
r <- summary(K0, show="dtb", save ="k")
print(r)
}
\seealso{
  Classes where this method applies:
  \code{\linkS4class{designkey}},
  \code{\linkS4class{keymatrix}},
  \code{\linkS4class{keyring}},
  \code{\linkS4class{listofdesignkeys}},
  \code{\linkS4class{listofkeyrings}},
  \code{\linkS4class{planordesign}}
}
\keyword{methods}
