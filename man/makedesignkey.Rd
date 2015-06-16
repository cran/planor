\name{makedesignkey}
\alias{makedesignkey}
\title{
  Turn Integer Matrices into an Object of Class designkey
}
\description{
  Create an object of class
  \code{\link[=designkey-class]{designkey}}  from a list of integer
  matrices.
}
\usage{
makedesignkey(keys, primes)
}
\arguments{
  \item{keys}{
a list of \code{n} integer matrices with column names.
}
  \item{primes}{
a vector of \code{n} prime numbers.
}
}
\details{
The names of the factors are extracted from the matrix column names.
}
\value{
An object of class \code{\link[=designkey-class]{designkey}}.
}
\author{Monod, H. }
\examples{
  mat1 <- cbind(diag(3),1)
  colnames(mat1) <- c("A","B","C","D")
  mat2 <- cbind(diag(2),c(1,2))
  colnames(mat2) <- c("E","F","G")
  mat.dk <- makedesignkey(list(mat1,mat2), primes=c(2,3))
  print(mat.dk)
  summary(mat.dk)
  alias(mat.dk)
  mat.plan <- planor.design(mat.dk)
}
