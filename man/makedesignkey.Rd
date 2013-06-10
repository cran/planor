\name{makedesignkey}
\alias{makedesignkey}
\title{
  A function to turn integer matrices into an object of class designkey
}
\description{
  Creates an object of class designkey directly from a list of integer
  matrices
}
\usage{
makedesignkey(keys, primes)
}
\arguments{
  \item{keys}{
a list of \code{n} integer matrices with column names
}
  \item{primes}{
a vector of \code{n} prime numbers
}
}
\details{
the names of the factors are extracted from the matrix column names
}
\value{
an object of class designkey
}
\author{
H. Monod
}
\seealso{
Class designkey
}
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
