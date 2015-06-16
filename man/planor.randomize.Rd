\name{planor.randomize}
\alias{planor.randomize}

\title{Randomize a Factorial Design from an Orthogonal Block Structure}

\usage{planor.randomize(blockformula, data, out.order, keep.initial=FALSE)}

\description{Randomize a factorial design according to a
  specified block structure formula.}

\arguments{\item{blockformula}{the block structure formula.}
  \item{data}{a data frame.}
  \item{out.order}{a list of \code{data} factors that will be used to
    order the rows of the randomized design; if missing, the factors of the
    block formula are used.}
  \item{keep.initial}{if TRUE, the initial row order of the design
    is stored in column \code{InitialUNITS} of the returned dataframe.}}

\value{The input data frame after randomization.}

\note{Each name in \code{blockformula} must correspond to a factor
  of the dataframe \code{data}. The only exception is
  \code{UNITS}. If \code{UNITS} is used in \code{blockformula}
  but absent from \code{data}, a factor is added to \code{data},
  with one level per row. See the examples below for the usage of
  \code{UNITS} in \code{blockformula}.}

\references{Bailey, R.A. (1983) Generalized wreath products of
  permutation groups. \emph{Proc. London Math. Soc.}, \bold{47}, 69--82.
    
  Kobilinsky A. (1989) Randomization of a cartesian block
  structure. Technical Report. Laboratoire de Biométrie. INRA.
  Versailles.}

\examples{
### Block design
Design <- data.frame(block=rep(1:4,rep(2,4)),
   treatment=c("A1","B1","A2","B2","A3","B3","A4","B4"))
### No within-block randomization
planor.randomize(~block, data=Design)
### Blocks and units within blocks randomization     
planor.randomize(~block/UNITS, data=Design) 
### Row-Column design
RowColDes <- data.frame(row=rep(1:3,rep(3,3)),col=rep(1:3,3),
treatment=LETTERS[c(1:3,2,3,1,3,1,2)],
oldRow=rep(1:3,rep(3,3)),oldCol=rep(1:3,3))
planor.randomize(~row*col, data=RowColDes)}

\author{Monod, H. and Bouvier, A. }

\keyword{design}
